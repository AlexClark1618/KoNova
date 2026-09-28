import numpy as np
import time
import os
import sys
import importlib
import matplotlib.pyplot as plt
from KoNova_Data_Filter import FilteredData

# ----------------------------- Initial Variables -----------------------------
HERE        = os.path.dirname(os.path.abspath(__file__))
MAP_FOLDER  = os.path.join(HERE, 'Bar_CH_Mappings')  # holds bar_ch_map_norm.csv, bar_ch_map_inverted.csv
BAR_SEP     = 16.5                       # mm, bar centre-to-centre
DELTA_Z     = 11.625 * 2.54 * 10          # mm, x-top to x-bottom (measured 4/29/2026)
RUN_SECONDS = 60
N_BARS      = 64
delta_x = 0#0.01 *10 #Dan's measured shifts between layers and stretch in bars
delta_y = 0#0.35 *10 
delta_x_stretch = 1#1.003 
delta_y_stretch = 1#0.9965 


FIDUCIAL_RADIUS = 1.0065e3 / 2          # mm # Drop edge bars
CENTER = ((N_BARS+1) / 2 * BAR_SEP,  # geometric centre of the bar array
          (N_BARS+1) / 2 * BAR_SEP)

# ----------------------------- Folder and File Handling -----------------------------
def folder_reader(folder_path, file_max = None):
    """Returns a list of filenames in folder_path matching run_name and ending with '_coinc.dat'."""
    coincidence_files = []
    for filename in os.listdir(folder_path):
        if filename.endswith('_coinc.dat'): #filename.startswith(run_name) and filename.endswith('_coinc.dat'):
            coincidence_files.append(os.path.join(folder_path, filename))
    if file_max:
        return coincidence_files[:file_max]
    else:
        return coincidence_files

def read_coincidence_file(coincidence_files):
    """Parse the .dat file into a list of event blocks: [[ts, qdc, ch], ...]."""
    events = []
    for file in coincidence_files:
        print(f"Processing file: {file}")

        with open(file, 'r') as f:
            while (header := f.readline()):
                n = sum(map(int, header.split()))
                try:
                    block = [
                        [int(ts), float(qdc), int(ch)]
                        for ts, qdc, ch in (f.readline().split() for _ in range(n))
                    ]
                    events.append(block)
                except ValueError:
                    print(f"Warning: Skipping malformed block in file {file}.")
                    continue

    return events

def load_detector_config(detector_name):
    """Import <detector_name>_initialization.py and return it (e.g. 'KN1' -> KN1_initialization)."""
    if HERE not in sys.path:
        sys.path.insert(0, HERE)
    config = importlib.import_module(f'{detector_name}_initialization')
    print(f"Loaded configuration for detector {config.DETECTOR_NAME}")
    for layer in sorted(config.OFFSETS):
        print(f"  Layer {layer}: offset {config.OFFSETS[layer]}, map {config.LAYER_MAPS[layer]}")
    return config

def read_bar_map(map_name, map_folder=MAP_FOLDER):
    """Read one bar/channel mapping CSV into {bar: base_ch}."""
    bar_map = {}
    with open(os.path.join(map_folder, map_name), 'r') as f:
        next(f)  # skip header
        for line in f:
            if not line.strip():
                continue
            bar, base_ch = map(int, line.strip().split(','))
            bar_map[bar] = base_ch
    return bar_map

def build_ch_to_bar(layer_maps, offsets, map_folder=MAP_FOLDER):
    """Map scaled channel number -> (layer, bar), using each layer's own mapping file."""
    ch_to_bar = {}
    cache = {}  # map filename -> {bar: base_ch}, so a shared map is only read once
    for layer, off in offsets.items():
        map_name = layer_maps[layer]
        if map_name not in cache:
            cache[map_name] = read_bar_map(map_name, map_folder)
        for bar, base_ch in cache[map_name].items():
            ch_to_bar[base_ch + off] = (layer, bar)
    return ch_to_bar

def write_to_stats_file():
    """Write the statistics to a text file."""
    stats_file_path = os.path.join(SAVE_FOLDER, f"{RUN_NAME}_stats.txt")
    with open(stats_file_path, 'w') as f:
        f.write(f"Run Name: {RUN_NAME}\n")
        #f.write(f"Total PETsys coincidence events: {len(events)}\n")
        #f.write(f"Raw PETsys Coincidence rate: {len(events) / (len(coincidence_files) * RUN_SECONDS)} Hz\n")
        #f.write(f"Filtered out {skipped} blocks with insufficient layers hit out of a total of {len(decoded)} blocks.\n")
        #f.write(f"Pruned {pruned_hits} stray bar hits from {pruned_blocks} blocks.\n")
        #f.write(f"Filtered out {too_wide} blocks with a cluster wider than {max_bars_per_layer} bars.\n")
        #f.write(f"Filtered out {all_single} blocks with only isolated single-bar hits in a layer.\n")
        #f.write(f"Filtered out {tied_multi} blocks with two or more multi-bar clusters in a layer.\n")
        #f.write(f"Blocks surviving adjacency cut: {len(kept_blocks)} of {len(filtered)}\n")
    print(f"Statistics written to {stats_file_path}")
# ------------------------ Decode to (layer, bar) ------------------
def decode_events(events, ch_to_bar):
    """Convert each [ts, qdc, ch] hit -> [layer, bar, qdc, ts], dropping unmapped channels."""
    decoded = []
    for block in events: #block = coincidence event
        decoded.append([
            [layer, bar, qdc, ts]
            for ts, qdc, ch in block
            if ch in ch_to_bar
            for (layer, bar) in [ch_to_bar[ch]]
        ])
    return decoded


'''
#From here I should cut any blocks without at least all 4 layers hit
def filter_decoded_events(decoded):
    filtered = []
    skipped = 0
    for block in decoded:
        layers_hit = {hit[0] for hit in block}
        if len(set(layers_hit)) == 4:
            #At least one hit in each layer
            filtered.append(block)
        else:
            skipped += 1

    print(f"Filtered out {skipped} blocks with insufficient layers hit out of a total of {len(decoded)} blocks.")
    return filtered

#With the filtered data there can still be hit events on a layer that are not adjacent.
#Rather than cutting the whole event, split the layer's bars into runs of consecutive bars
#(e.g. 3,4,57 -> [3,4] and [57]), keep the run that looks like the real hit, and drop the rest as noise.

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


#Reason for picking a cluster by size alone: two neighbouring bars firing together by chance is
#very unlikely, so a multi-bar cluster is a real hit while a lone bar off on its own is consistent
#with noise. QDC is not used anywhere in this cut.

AMBIGUOUS_ALL_SINGLE = 'all_single'   #every cluster is one bar wide - nothing distinguishes them
AMBIGUOUS_TIED_MULTI = 'tied_multi'   #two or more multi-bar clusters - looks like two real hits


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

def bar_adjacency_cut(filtered, max_bars_per_layer=5):
    """Prune stray non-adjacent bar hits, then cut events with too wide a cluster."""
    kept_blocks = []
    too_wide = 0
    all_single = 0
    tied_multi = 0
    pruned_blocks = 0
    pruned_hits = 0

    for block in filtered:
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

    print(f"Pruned {pruned_hits} stray bar hits from {pruned_blocks} blocks.")
    print(f"Filtered out {too_wide} blocks with a cluster wider than {max_bars_per_layer} bars.")
    print(f"Filtered out {all_single} blocks with only isolated single-bar hits in a layer.")
    print(f"Filtered out {tied_multi} blocks with two or more multi-bar clusters in a layer.")
    print(f"Blocks surviving adjacency cut: {len(kept_blocks)} of {len(filtered)}")
    return kept_blocks

def filter_by_time(kept_blocks, max_time_ns=10):

    qdc_fail_filtered_hits = []
    qdc_pass_filtered_hits = []

    time_filtered_blocks = []

    for block in kept_blocks:
        temp_block =[]
        ts_list = [hit[3] for hit in block]
        if (max(ts_list) - min(ts_list))/1000 > max_time_ns:
            # print(f'Ts spread greater than {max_time_ns}ns')
            for hit in block:
                if (hit[3] - min(ts_list))/1000 <= max_time_ns:
                    temp_block.append(hit)
                    qdc_pass_filtered_hits.append(hit[2])

                else:
                    qdc_fail_filtered_hits.append(hit[2])
            #print("old:",block)
            #print("new:",temp_block)

            time_filtered_blocks.append(temp_block)

        else:
            time_filtered_blocks.append(block)

    plt.hist(qdc_fail_filtered_hits, bins='fd', color='red', alpha=0.5, label=f"Time Cut Failures\nMedian: {np.median(qdc_fail_filtered_hits):.2f} | Std: {np.std(qdc_fail_filtered_hits):.2f}")
    plt.hist(qdc_pass_filtered_hits, bins='fd', color='blue', alpha=0.5, label=f"Time Cut Passes\nMedian: {np.median(qdc_pass_filtered_hits):.2f} | Std: {np.std(qdc_pass_filtered_hits):.2f}")

    plt.title(f'QDC Distribution for Time Cut Passes and Failures | Max Time Spread: {max_time_ns} ns')
    plt.xlabel('QDC Value')
    plt.ylabel('Frequency')
 
    plt.legend()
    plt.show()
    return time_filtered_blocks
'''

def per_layer_stats(kept_blocks):
    """Return (freq, qdc) dicts keyed by layer."""
    freq_dist = {l: {b: 0 for b in range(1, N_BARS + 1)} for l in range(1, 5)}
    qdc_dist  = {l: [] for l in range(1, 5)}
    time_dist  = {l: [] for l in range(1, 5)}

    ts_greater_than_15ns = 0
    for block in kept_blocks:
        ts_list = [hit[3] for hit in block]
        if (max(ts_list) - min(ts_list))/1000 > 20:
            ts_greater_than_15ns += 1
            #print(f'Ts spread greater than 15ns {(max(ts_list) - min(ts_list))/1000} ns')
            #print(sorted(ts_list))
        for layer, bar, qdc_value, _ in block:

            freq_dist[layer][bar] += 1
            qdc_dist[layer].append(qdc_value)
            time_dist[layer].append((max(ts_list) - min(ts_list))/1000)


    print(f"Blocks with time spread greater than 15ns: {ts_greater_than_15ns} of {len(kept_blocks)}")
    return freq_dist, qdc_dist, time_dist

def bar_frequency_and_qdc_distribution_plots(data, graph):
    freq, qdc, ts = per_layer_stats(data)

    for layer in range(1, 5):
        plt.bar(freq[layer].keys(), freq[layer].values())
        mean_freq = np.mean(list(freq[layer].values()))
        std_freq = np.std(list(freq[layer].values()))
        plt.title(f'Layer {layer} Bar Hit Frequency\nMean: {mean_freq:.2f} | Std: {std_freq:.2f}')
        plt.xlabel('Bar Number')
        plt.ylabel('Frequency')

        filename = f"{RUN_NAME}_Layer_{layer}_Bar_Hit_Frequency"
        print(f'{filename} saved')
        filepath = os.path.join(SAVE_FOLDER, filename)
        plt.savefig(filepath, dpi=300)

        if graph:
            plt.show()
        
        plt.clf()   

        plt.hist(qdc[layer], bins='fd')
        mean_qdc = np.mean(qdc[layer])
        std_qdc = np.std(qdc[layer])
        plt.title(f'Layer {layer} QDC Distribution\nMean: {mean_qdc:.2f} | Std: {std_qdc:.2f}')
        plt.xlim(-1, 10)
        plt.xlabel('QDC Value')
        plt.ylabel('Frequency')

        filename = f"{RUN_NAME}_Layer_{layer}_QDC_Distribution"
        print(f'{filename} saved')
        filepath = os.path.join(SAVE_FOLDER, filename)
        plt.savefig(filepath, dpi=300)

        if graph:
            plt.show()
        
        plt.clf() 

        plt.hist(ts[layer], bins='fd')
        mean_ts = np.mean(ts[layer])
        std_ts = np.std(ts[layer])
        plt.title(f'Layer {layer} Time Distribution\nMean: {mean_ts:.2f} | Std: {std_ts:.2f}')
        plt.xlabel('Time (ns)')
        plt.xlim(0, 100)
        plt.ylabel('Frequency')
        plt.yscale('log')

        filename = f"{RUN_NAME}_Layer_{layer}_Time_Distribution"
        print(f'{filename} saved')
        filepath = os.path.join(SAVE_FOLDER, filename)
        plt.savefig(filepath, dpi=300)

        if graph:
            plt.show()
        
        plt.clf() 
#--- Original whole-block adjacency cut, kept for reference ---------------------------
#Cuts the entire event whenever a layer has >3 bars or any non-adjacent pair, instead of
#pruning the stray hit. Note: it mutates `filtered` while iterating over it, so the loop
#skips elements and the counts undercount - re-check the numbers before trusting them.
#Needs `import numpy as np`.
#
# def bar_adjacency_cut(filtered):
#     greater_than_3_bars = 0
#     greater_than_1_non_adjacent = 0
#     at_least_1_adjacent = 0
#     for block in filtered:
#         by_layer = {l: [h[1] for h in block if h[0] == l] for l in range(1, 5)}
#
#         for layer, bars in by_layer.items():
#             #I think the max bars hit per layer should be set to 3. If two are adjacent and one is not it should be allowed, but with the none adjacent bar cut
#             #print(f"Layer {layer} bars: {bars}")
#             if len(bars)>3:
#                 #print(f"Warning: More than 3 bars hit in layer {layer}. Bars: {sorted(bars)}")
#
#                 try:
#                     filtered.remove(block)
#                     greater_than_3_bars += 1
#
#                 except ValueError:
#                     pass
#
#             elif len(bars)>1 and np.diff(sorted(bars)).max() > 1:
#
#                 #print(f"Warning: Non-adjacent bars hit in layer {layer}. Bars: {sorted(bars)}")
#
#                 try:
#                     filtered.remove(block)
#                     greater_than_1_non_adjacent += 1
#                     if np.diff(sorted(bars)).min() == 1:
#                         at_least_1_adjacent += 1
#
#                 except ValueError:
#                     pass
#
#     print("Filtered filteres", len(filtered))
#     print(f"Filtered out {greater_than_3_bars} blocks with more than 3 bars hit in a layer.")
#     print(f"Filtered out {greater_than_1_non_adjacent} blocks with non-adjacent bars hit in a layer.")
#     print(f"Blocks with at least 1 adjacent bar: {at_least_1_adjacent}")
#--------------------------------------------------------------------------------------


def main(config):

    coincidence_files = folder_reader(SUB_DATA_FOLDER_PATH, file_max=1)
    print(f"Found {len(coincidence_files)} files for run {RUN_NAME}.")
    events    = read_coincidence_file(coincidence_files)

    print(f"Total PETsys coincidence events: {len(events)}")
    print(f"Raw PETsys Coincidence rate: {len(events) / (len(coincidence_files) * RUN_SECONDS)} Hz")

    ch_to_bar = build_ch_to_bar(config.LAYER_MAPS, config.OFFSETS)
    raw_data   = decode_events(events, ch_to_bar)
    filtered_data = FilteredData(raw_data,
                 max_bars_per_layer=5,
                 max_time_ns=100,
                 plot_time_cut=True,
                 save_folder=SAVE_FOLDER,
                 run_name=RUN_NAME,
                 show_plots=False)


    bar_frequency_and_qdc_distribution_plots(filtered_data, graph= False)

    write_to_stats_file()

if __name__ == '__main__':
    #FILE_PATH   = 'KNVA-20260514-01-00079_coinc.dat'

    DETECTOR    = 'KN2'   # 'KN1' or 'KN2' -> loads <DETECTOR>_initialization.py
    config      = load_detector_config(DETECTOR)

    DATA_FOLDER_PATH = r"C:\\Users\\aclark2\\Desktop\\KoNova\\PETsys Data"
    SAVE_FOLDER_PATH = r"C:\\Users\\aclark2\\Desktop\\KoNova\\PETsys Plots"
    RUN_NAME    = 'KN2_Lab_Test'
    SAVE_RUN_NAME = 'KN2_Lab_Test' #RUN_NAME
    SUB_DATA_FOLDER_PATH = os.path.join(DATA_FOLDER_PATH, RUN_NAME)
    SAVE_FOLDER = os.path.join(SAVE_FOLDER_PATH, SAVE_RUN_NAME)
    os.makedirs(SAVE_FOLDER, exist_ok=True) 
    print(f'Save directory created {SAVE_FOLDER}')

    main(config)