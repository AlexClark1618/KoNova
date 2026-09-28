import os
import numpy as np
import matplotlib.pyplot as plt

#Reason for picking a cluster by size alone: two neighbouring bars firing together by chance is
#very unlikely, so a multi-bar cluster is a real hit while a lone bar off on its own is consistent
#with noise. QDC is not used anywhere in this cut.

AMBIGUOUS_ALL_SINGLE = 'all_single'   #every cluster is one bar wide - nothing distinguishes them
AMBIGUOUS_TIED_MULTI = 'tied_multi'   #two or more multi-bar clusters - looks like two real hits


def cluster_bars(bars):
    """Split bar numbers into clusters of consecutive bars. [3,4,57] -> [[3,4],[57]]."""
    clusters = []
    for bar in sorted(bars):
        if clusters and bar - clusters[-1][-1] <= 1:   # <=1 so a repeated bar joins its own cluster
            if bar != clusters[-1][-1]:
                clusters[-1].append(bar)
        else:
            clusters.append([bar])
    return clusters


def prune_layer(hits):
    """Keep only the widest cluster of consecutive bars in one layer's hits.

    Returns (kept_hits, n_dropped, ambiguity). `ambiguity` is None when one cluster wins
    outright, otherwise a reason string saying why the layer could not be resolved.
    """
    clusters = cluster_bars(h[1] for h in hits)
    if len(clusters) == 1:
        return hits, 0, None

    widest = max(len(c) for c in clusters)
    winners = [c for c in clusters if len(c) == widest]

    if len(winners) > 1:
        return hits, 0, AMBIGUOUS_ALL_SINGLE if widest == 1 else AMBIGUOUS_TIED_MULTI

    #A single one-bar cluster can only win if it is the sole cluster, which is handled above.
    members = set(winners[0])
    kept = [h for h in hits if h[1] in members]
    return kept, len(hits) - len(kept), None


class FilteredData:
    """Run the full cut chain on decoded events.

    Usage:
        from KoNova_Data_Filter import FilteredData

        cuts = FilteredData(decoded, max_bars_per_layer=5, max_time_ns=15)
        clean = cuts.blocks          # blocks surviving every cut

    Intermediate stages stay available as .layer_filtered / .adjacency_filtered /
    .time_filtered, and per-cut counts in .stats.
    """

    def __init__(self, decoded, max_bars_per_layer=5, max_time_ns=10,
                 verbose=True, plot_time_cut=False, save_folder=None,
                 run_name='', show_plots=False, run=True):
        self.decoded = decoded
        self.max_bars_per_layer = max_bars_per_layer
        self.max_time_ns = max_time_ns
        self.verbose = verbose
        self.plot_time_cut = plot_time_cut
        self.save_folder = save_folder
        self.run_name = run_name
        self.show_plots = show_plots

        self.layer_filtered = None
        self.adjacency_filtered = None
        self.time_filtered = None
        self.stats = {}

        if run:
            self.run()

    # ------------------------------------------------------------------
    @property
    def blocks(self):
        """Blocks left after every cut (falls back to the last stage that ran)."""
        for stage in (self.time_filtered, self.adjacency_filtered,
                      self.layer_filtered, self.decoded):
            if stage is not None:
                return stage

    def __len__(self):
        return len(self.blocks)

    def __iter__(self):
        return iter(self.blocks)

    def _log(self, message):
        if self.verbose:
            print(message)

    # ------------------------------------------------------------------
    def run(self):
        """Apply every cut in order and return the surviving blocks."""
        self.filter_by_layer()
        self.filter_by_bar_adjacency()
        self.filter_by_time()
        return self.blocks

    def filter_by_layer(self):
        """Filter out blocks that do not have at least one hit in each of the four layers."""
        filtered = []
        skipped = 0
        for block in self.decoded:
            layers_hit = {hit[0] for hit in block}
            if len(layers_hit) == 4:
                filtered.append(block)
            else:
                skipped += 1

        self._log(f"Filtered out {skipped} blocks with insufficient layers hit "
                  f"out of a total of {len(self.decoded)} blocks.")

        self.stats['layer_skipped'] = skipped
        self.stats['layer_kept'] = len(filtered)
        self.layer_filtered = filtered
        return filtered

    def filter_by_bar_adjacency(self):
        """Prune stray non-adjacent bar hits, then cut events with too wide a cluster.

        With the filtered data there can still be hit events on a layer that are not
        adjacent. Rather than cutting the whole event, split the layer's bars into runs
        of consecutive bars (e.g. 3,4,57 -> [3,4] and [57]), keep the run that looks
        like the real hit, and drop the rest as noise.
        """
        if self.layer_filtered is None:
            self.filter_by_layer()

        max_bars_per_layer = self.max_bars_per_layer
        kept_blocks = []
        too_wide = 0
        all_single = 0
        tied_multi = 0
        pruned_blocks = 0
        pruned_hits = 0

        for block in self.layer_filtered:
            new_block = []
            dropped_here = 0
            reject = False

            for layer in range(1, 5):
                hits = [h for h in block if h[0] == layer]
                if not hits:
                    continue

                kept, n_dropped, ambiguity = prune_layer(hits)

                if ambiguity is not None:
                    if ambiguity == AMBIGUOUS_ALL_SINGLE:
                        all_single += 1
                    else:
                        tied_multi += 1
                    reject = True
                    break

                #After pruning, the surviving bars are consecutive by construction, so this is
                #purely a cluster-width cut: a shower spread over too many bars isn't a clean track.
                if len({h[1] for h in kept}) > max_bars_per_layer:
                    too_wide += 1
                    reject = True
                    break

                dropped_here += n_dropped
                new_block.extend(kept)

            if reject:
                continue

            if dropped_here:
                pruned_blocks += 1
                pruned_hits += dropped_here
            kept_blocks.append(new_block)

        self._log(f"Pruned {pruned_hits} stray bar hits from {pruned_blocks} blocks.")
        self._log(f"Filtered out {too_wide} blocks with a cluster wider than {max_bars_per_layer} bars.")
        self._log(f"Filtered out {all_single} blocks with only isolated single-bar hits in a layer.")
        self._log(f"Filtered out {tied_multi} blocks with two or more multi-bar clusters in a layer.")
        self._log(f"Blocks surviving adjacency cut: {len(kept_blocks)} of {len(self.layer_filtered)}")

        self.stats.update({
            'pruned_hits': pruned_hits,
            'pruned_blocks': pruned_blocks,
            'adjacency_too_wide': too_wide,
            'adjacency_all_single': all_single,
            'adjacency_tied_multi': tied_multi,
            'adjacency_kept': len(kept_blocks),
        })
        self.adjacency_filtered = kept_blocks
        return kept_blocks

    def filter_by_time(self):
        """Drop hits whose timestamp sits more than max_time_ns after the block's first hit."""
        if self.adjacency_filtered is None:
            self.filter_by_bar_adjacency()

        max_time_ns = self.max_time_ns
        qdc_fail_filtered_hits = []
        qdc_pass_filtered_hits = []
        time_filtered_blocks = []
        blocks_trimmed = 0

        for block in self.adjacency_filtered:
            ts_list = [hit[3] for hit in block]
            if (max(ts_list) - min(ts_list)) / 1000 > max_time_ns:
                blocks_trimmed += 1
                temp_block = []
                t0 = min(ts_list)
                for hit in block:
                    if (hit[3] - t0) / 1000 <= max_time_ns:
                        temp_block.append(hit)
                        qdc_pass_filtered_hits.append(hit[2])
                    else:
                        qdc_fail_filtered_hits.append(hit[2])

                time_filtered_blocks.append(temp_block)
            else:
                time_filtered_blocks.append(block)

        self._log(f"Trimmed hits from {blocks_trimmed} blocks with a time spread over {max_time_ns} ns "
                  f"({len(qdc_fail_filtered_hits)} hits dropped, {len(qdc_pass_filtered_hits)} kept).")

        self.stats.update({
            'time_blocks_trimmed': blocks_trimmed,
            'time_hits_dropped': len(qdc_fail_filtered_hits),
            'time_kept': len(time_filtered_blocks),
        })
        self.time_filtered = time_filtered_blocks

        if self.plot_time_cut:
            self._plot_time_cut_qdc(qdc_fail_filtered_hits, qdc_pass_filtered_hits)

        return time_filtered_blocks

    # ------------------------------------------------------------------
    def _plot_time_cut_qdc(self, fail_hits, pass_hits):
        """QDC distribution of the hits the time cut dropped vs. kept."""
        if not fail_hits and not pass_hits:
            self._log("No blocks exceeded the time spread - skipping the time cut QDC plot.")
            return

        if fail_hits:
            plt.hist(fail_hits, bins='fd', color='red', alpha=0.5,
                     label=f"Time Cut Failures\nMedian: {np.median(fail_hits):.2f} | Std: {np.std(fail_hits):.2f}")
        if pass_hits:
            plt.hist(pass_hits, bins='fd', color='blue', alpha=0.5,
                     label=f"Time Cut Passes\nMedian: {np.median(pass_hits):.2f} | Std: {np.std(pass_hits):.2f}")

        plt.title(f'QDC Distribution for Time Cut Passes and Failures | Max Time Spread: {self.max_time_ns} ns')
        plt.xlabel('QDC Value')
        plt.ylabel('Frequency')
        plt.legend()

        if self.save_folder:
            filename = f"{self.run_name}_Time_Cut_QDC_Distribution.png"
            filepath = os.path.join(self.save_folder, filename)
            plt.savefig(filepath, dpi=300)
            self._log(f'{filename} saved')

        if self.show_plots:
            plt.show()

        plt.clf()

    # ------------------------------------------------------------------
    def summary(self):
        """One-line-per-stage rundown of how many blocks survived each cut."""
        lines = [f"Decoded blocks:            {len(self.decoded)}"]
        if self.layer_filtered is not None:
            lines.append(f"After 4-layer cut:         {len(self.layer_filtered)}")
        if self.adjacency_filtered is not None:
            lines.append(f"After bar adjacency cut:   {len(self.adjacency_filtered)}")
        if self.time_filtered is not None:
            lines.append(f"After time cut:            {len(self.time_filtered)}")
        return "\n".join(lines)
