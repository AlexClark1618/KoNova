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
RUN_SECONDS = 600
N_BARS      = 64
delta_x = 0#0.01 *10 #Dan's measured shifts between layers and stretch in bars
delta_y = 0#0.35 *10 
delta_x_stretch = 1#1.003 
delta_y_stretch = 1#0.9965 


FIDUCIAL_RADIUS = 1.0395e3 / 2          # mm # Drop edge bars
CENTER = ((N_BARS+1) / 2 * BAR_SEP,  # geometric centre of the bar array
          (N_BARS+1) / 2 * BAR_SEP)

# ----------------------------- Folder and File Handling -----------------------------
def folder_reader(folder_path, file_max = None, file_size_in_MB=None):
    """Returns a list of filenames in folder_path matching run_name and ending with '_coinc.dat'."""
    coincidence_files = []
    for filename in os.listdir(folder_path):
        if filename.endswith('_coinc.dat'): #filename.startswith(run_name) and filename.endswith('_coinc.dat'):
            coincidence_files.append(os.path.join(folder_path, filename))
    if file_max:
        return coincidence_files[:file_max]
    if file_size_in_MB:
        for f in coincidence_files:
            if os.path.getsize(f) < file_size_in_MB * 1000 * 1000 or file_size_in_MB :
                print(f"Warning: File {f} is smaller than {file_size_in_MB} MB and will be skipped.")
        return [f for f in coincidence_files if os.path.getsize(f) > file_size_in_MB * 1000 * 1000]
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

def write_to_stats_file(cuts, n_events, n_files, tracks_full, tracks_fiducial):
    """Write the run's event counts, per-cut tallies and rates to a text file."""
    total_seconds = n_files * RUN_SECONDS
    s = cuts.stats
    stats_file_path = os.path.join(SAVE_FOLDER, f"{RUN_NAME}_stats.txt")
    with open(stats_file_path, 'w') as f:
        f.write(f"Run Name: {RUN_NAME}\n")
        f.write(f"Files processed: {n_files} ({total_seconds} s)\n")
        f.write(f"Total PETsys coincidence events: {n_events}\n")
        f.write(f"Raw PETsys Coincidence rate: {n_events / total_seconds} Hz\n")
        f.write(f"Filtered out {s['layer_skipped']} blocks with insufficient layers hit out of a total of {len(cuts.decoded)} blocks.\n")
        f.write(f"Pruned {s['pruned_hits']} stray bar hits from {s['pruned_blocks']} blocks.\n")
        f.write(f"Filtered out {s['adjacency_too_wide']} blocks with a cluster wider than {cuts.max_bars_per_layer} bars.\n")
        f.write(f"Filtered out {s['adjacency_all_single']} blocks with only isolated single-bar hits in a layer.\n")
        f.write(f"Filtered out {s['adjacency_tied_multi']} blocks with two or more multi-bar clusters in a layer.\n")
        f.write(f"Blocks surviving adjacency cut: {s['adjacency_kept']} of {s['layer_kept']}\n")
        f.write(f"Trimmed hits from {s['time_blocks_trimmed']} blocks with a time spread over {cuts.max_time_ns} ns "
                f"({s['time_hits_dropped']} hits dropped).\n")
        f.write(f"Tracks reconstructed (full area): {len(tracks_full)}\n")
        f.write(f"Tracks reconstructed (fiducial): {len(tracks_fiducial)}\n")
        f.write(f"Full Area Coincidence Rate: {len(tracks_full) / total_seconds} Hz\n")
        f.write(f"Fiducial Area Coincidence Rate: {len(tracks_fiducial) / total_seconds} Hz\n")
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


# --------------------------- Diagnostics --------------------------
def per_layer_stats(kept_blocks):
    """Return (freq, qdc) dicts keyed by layer."""
    freq_dist = {l: {b: 0 for b in range(1, N_BARS + 1)} for l in range(1, 5)}
    qdc_dist  = {l: [] for l in range(1, 5)}
    time_dist  = {l: [] for l in range(1, 5)}

    ts_greater_than_20ns = 0
    for block in kept_blocks:
        ts_list = [hit[3] for hit in block]
        if (max(ts_list) - min(ts_list))/1000 > 20:
            ts_greater_than_20ns += 1

        for layer, bar, qdc_value, _ in block:

            freq_dist[layer][bar] += 1
            qdc_dist[layer].append(qdc_value)
            time_dist[layer].append((max(ts_list) - min(ts_list))/1000)

    print(f"Blocks with time spread greater than 20ns: {ts_greater_than_20ns} of {len(kept_blocks)}")
    return freq_dist, qdc_dist, time_dist

from collections import Counter

def bar_hits_per_layer(decoded, cuts=True, graph=False):
    num_bar_hit_freq_per_layer = {l: [] for l in range(1, 5)}

    if cuts:
        name_add_on = "w Cuts"
    else:
        name_add_on = "w-o Cuts"

    summary_filepath = os.path.join(SAVE_FOLDER, f"{RUN_NAME}_Bar_Hit_Multiplicity_Summary_({name_add_on}).txt")

    for block in decoded:
        
        if cuts:
            by_layer = {l: [h for h in block if h[0] == l] for l in range(1, 5)}

            if any(len(by_layer[l]) == 0 for l in range(1, 5)):
                continue

            max_diff_flag = False
            for l in range(1, 5):
                bars = np.array(sorted(h[1] for h in by_layer[l]))
                if len(bars) == 1:
                    continue
                if np.max(np.diff(bars)) > 1:
                    max_diff_flag = True
                    break

        if cuts and max_diff_flag:
            continue

        block = np.array(block)
        int_array = [int(x) for x in block[:, 0]]
        layer_hit_counts_per_event = Counter(int_array)

        # for layer, counts in layer_hit_counts_per_event.items():
        #     num_bar_hit_freq_per_layer[layer].append(counts)
        for layer in range(1, 5):
            counts = layer_hit_counts_per_event.get(layer, 0)  # 0 if layer not in counter
            num_bar_hit_freq_per_layer[layer].append(counts)
            
    with open(summary_filepath, 'w') as f:
        for layer, count_list in num_bar_hit_freq_per_layer.items():
            count_hist = Counter(count_list)
            #print(f"Layer {layer}: {count_hist}")
            f.write(f"Layer {layer}: {dict(sorted(count_hist.items()))}\n")

            plt.bar(count_hist.keys(), count_hist.values())
            plt.xlim(-0.5, 5.5)
            plt.title(f'Layer {layer} Bar Hit Multiplicity Distribution ({name_add_on})')

            filename = f"{RUN_NAME}_Layer_{layer}_Bar_Hit_Multiplicity_Distribution_({name_add_on}).png"
            filepath = os.path.join(SAVE_FOLDER, filename)
            print(f'{filename} saved')
            plt.savefig(filepath, dpi=300)

            if graph:
                plt.show()

            plt.clf()
# --------------------- Track / position building -----------------
def layer_position(hits, sep=BAR_SEP):
    """Mean bar position (mm) for a list of hits in one layer, or None if empty."""
    if not hits: 
        #If layer empty return None, should be caught by build_tracks and skipped
        return None
    bars = np.array([h[1] for h in hits])

    deadzone = 0 #mm (estimate)
    if len(set(bars)) > 1: 
        #If multiple bars hit, return random position between lowest and highest hit bars hit
        low  = float((bars.min() * sep) - BAR_SEP + deadzone)   # Center of lowest hit bar
        high = float((bars.max() * sep) + BAR_SEP - deadzone)  # Center of highest hit bar
        bar_pos = np.random.uniform(low, high)
        return bar_pos
    
    else:
        low  = float((bars * sep) - BAR_SEP)   # Low of single bar hit with deadzone
        high = float((bars * sep) + BAR_SEP)  # High of single bar hit with deadzone
        bar_pos = np.random.uniform(low, high)
        return bar_pos
        #return int(bars) * sep # Bar 1 center = 16.5, Between Bar 1 and 2 = 24.75, ... Bar 64 center = 1056.0


def bars_adjacent(layer_hits):
    """True if 1 hit, or if all hit bars in the layer are consecutive."""
    if len(layer_hits) == 1:
        #print(f"Only one hit in layer, accepting by default {layer_hits}")
        return True
    bars = np.array(sorted(h[1] for h in layer_hits))
    #print(f"Checking adjacency for bars: {bars}")
    #print(f"Bar differences: {np.diff(bars)}")

    if np.max(np.diff(bars)) > 1:
        pass
        #print(f"Non-adjacent bars found: {bars}")

    #print(bars[0], bars[-1], len(bars))
    
    return bars[-1] - bars[0] == len(bars) - 1 and len(set(bars)) == len(bars)

def build_tracks(decoded, accept=bars_adjacent, cut_edge_bars=True):
    """
    Build [(x_top, y_top), (x_bottom, y_bottom)] per event.

    accept : optional callable(layer_hits) -> bool, applied per layer.
             Event is dropped if any required layer fails or is empty.
    """
    tracks = []
    dxdy=[]
    skipped_none = 0
    skipped_adjacency = 0
    skipped_edge = 0
    for block in decoded:
        by_layer = {l: [h for h in block if h[0] == l] for l in range(1, 5)}
        #print(by_layer)
        if any(len(by_layer[l]) == 0 for l in range(1, 5)): #Skip if any layer is empty
            #print(f"Empty layer , skipping event.")
            #time.sleep(3)
            skipped_none += 1
            continue
             
        if cut_edge_bars and any(
            any(h[1] == 1 or h[1] == N_BARS for h in by_layer[l])
            for l in range(1, 5)
        ):
            skipped_edge += 1
            continue

        elif accept and not all(accept(by_layer[l]) for l in range(1, 5)):
            #print(f"Non-adjacent bars in layer, skipping event. Layer hits: {by_layer}")
            skipped_adjacency += 1
            continue

        pos = {l: layer_position(by_layer[l]) for l in range(1, 5)}
        if any(p is None for p in pos.values()):
            #print(f"Could not determine position for event, skipping. Layer positions: {pos}")
            continue

        # layer 1=x_bottom, 2=y_bottom, 3=x_top, 4=y_top
        top    = (pos[3], pos[4])
        bottom = (pos[1], pos[2])
        
        #dx = bottom[0] - top[0]
        #dy = bottom[1] - top[1]
        
        #dxdy.append([dx,dy])
        
        tracks.append([top, bottom])
    print(f"Events skipped (empty layers): {skipped_none}")
    print(f"Events skipped (non-adjacent bars): {skipped_adjacency}")
    print(f"Events skipped (edge bars): {skipped_edge}")
    print(f"Events kept: {len(tracks)}")
    
    #dxdy= np.array(dxdy)
    
    #dx_hist, dy_hist = np.histogram(dxdy[:,0], bins=100), np.histogram(dxdy[:,1], bins=100)
    
    #plt.hist(dxdy[:,0], bins=100)
    #plt.show()
    #plt.hist(dxdy[:,1], bins=100)
    #plt.show()
    
    
    return tracks


# ----------------------------- Angles -----------------------------

both_zero = 0 
dx_zero = 0
dy_zero = 0
call_count = 0
def track_angles(p_top, p_bottom, delta_z=DELTA_Z):
    global both_zero, dx_zero, dy_zero, call_count
    """Return (zenith, azimuth) in degrees. Zenith from +z, azimuth CCW from +x."""
    dx = p_bottom[0] - p_top[0]
    dy = p_bottom[1] - p_top[1]
    call_count += 1
    
    if dx == 0 and dy == 0:
        both_zero += 1
    elif dx == 0:
        dx_zero += 1
    elif dy == 0:
        dy_zero += 1
    #print(f'dx: {dx}')
    v  = np.array([(dx * delta_y_stretch) + delta_y , (dy * delta_x_stretch) + delta_x, delta_z])
    #print(f"v[x]: {v[0]}")
    zenith  = np.degrees(np.arccos(np.clip(v[2] / np.linalg.norm(v), -1.0, 1.0)))
    azimuth = np.degrees(np.arctan2(v[1], v[0])) % 360
    
    return zenith, azimuth

def in_fiducial(point, center=CENTER, radius=FIDUCIAL_RADIUS):
    return (point[0] - center[0])**2 + (point[1] - center[1])**2 <= radius**2

def fiducial_filter(tracks):
    """Keep only tracks whose both endpoints lie within the cylinder."""
    return [t for t in tracks if in_fiducial(t[0]) and in_fiducial(t[1])]

def compute_angle_distributions(tracks, delta_z=DELTA_Z):
    if not tracks:
        raise ValueError(
            "No tracks survived build_tracks - check the skip counts above. "
            "If everything was skipped for empty layers, either the run was taken in "
            "2-fold coincidence (no event has all 4 layers, so no track can be built), "
            "or the DETECTOR config does not match the data (wrong channel offsets or "
            "bar/channel map)."
        )
    angles = np.array([track_angles(t, b, delta_z) for t, b in tracks])
    return angles[:, 0], angles[:, 1]   # zenith, azimuth

# ----------------------------- Plots -----------------------------

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

def zenith_and_azimuth_distribution_plot(zenith, azimuth, graph, full_area):

    if full_area:
        name_add_on = "Full Area"

    else:
        name_add_on = "Fiducial"

    #Zenith Plotting
    plt.hist(zenith, bins=90)
    mean_zenith = np.mean(zenith)
    std_zenith = np.std(zenith)
    plt.title(f'Zenith Angle Distribution {name_add_on}\nMean: {mean_zenith:.2f}° | Std: {std_zenith:.2f}°')
    plt.xlabel('Zenith (°)') 
    plt.ylabel('Counts')
    filename = f"{RUN_NAME}_Zenith_Angle_Distribution_{name_add_on}"
    print(f'{filename} saved')
    filepath = os.path.join(SAVE_FOLDER, filename)
    plt.savefig(filepath, dpi=300)

    if graph:
        plt.show()
    
    plt.clf() 

    #Asimuth Plotting
    plt.hist(azimuth, bins=90)
    plt.title(f'Azimuth Angle Distribution {name_add_on}')
    plt.xlabel('Azimuth (°)')
    plt.ylabel('Counts')
    filename = f"{RUN_NAME}_Azimuth_Angle_Distribution_{name_add_on}"
    print(f'{filename} saved')
    filepath = os.path.join(SAVE_FOLDER, filename)
    plt.savefig(filepath, dpi=300)

    if graph:
        plt.show()
    
    plt.clf() 

def anglular_heatmap(zenith, azimuth, graph, full_area):

    if full_area:
        name_add_on = "Full Area"

    else:
        name_add_on = "Fiducial"

    azimuth_to_radians = np.radians(azimuth)
    azimuth_bins = np.linspace(0, 2*np.pi, 65)  # 64 bins + 1 edge
    zenith_bins = np.linspace(0, 90, 10)        # 9 bins + 1 edge

    X1, X2 = np.meshgrid(azimuth_bins, zenith_bins)

    angular_hist = np.histogram2d(azimuth_to_radians, zenith, bins=[azimuth_bins, zenith_bins])

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.set_theta_zero_location('E')
    ax.set_theta_direction(1)
    mesh = ax.pcolormesh(X1.T, X2.T, angular_hist[0], cmap='viridis', shading='flat')
    fig.colorbar(mesh, ax=ax, pad=0.1,label='Counts')
    ax.set_title(f'Angular Distribution Heatmap {name_add_on}')

    filename = f"{RUN_NAME}_Angular_Distribution_Heatmap_{name_add_on}"
    print(f'{filename} saved')
    filepath = os.path.join(SAVE_FOLDER, filename)
    plt.savefig(filepath, dpi=300)

    if graph:
        plt.show()
    
    plt.clf() 

    #Solid angle correction:
    azimuth_delta = 2 * np.pi / (len(azimuth_bins)-1)
    print(f"Azimuth Delta: {azimuth_delta}")  
    
    solid_angle_correction_array = [azimuth_delta * (np.cos(np.radians(zenith_bins[i])) - np.cos(np.radians(zenith_bins[i+1]))) for i in range(len(zenith_bins)-1)]
    #print(f"Solid Angle: {solid_angle_correction_array}")

    flux_per_solid_angle = angular_hist[0] / np.array(solid_angle_correction_array)
    #print(f"Flux Density Per Steradian: {flux_per_solid_angle}")

    #Normalize by max to give relative intensity    
    flux_per_solid_angle /= np.max(flux_per_solid_angle)

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.set_theta_zero_location('E')
    ax.set_theta_direction(1)
    mesh = ax.pcolormesh(X1.T, X2.T, flux_per_solid_angle, cmap='viridis', shading='flat')
    fig.colorbar(mesh, ax=ax, pad=0.1, label='Relative Intensity (Normalized by Max)')
    ax.set_title(f'Angular Distribution Heatmap (Solid Angle Corrected) {name_add_on}')

    filename = f"{RUN_NAME}_Angular_Distribution_Heatmap_Solid_Angle_Corrected_{name_add_on}"
    print(f'{filename} saved')
    filepath = os.path.join(SAVE_FOLDER, filename)
    plt.savefig(filepath, dpi=300)

    if graph:
        plt.show()
    
    plt.clf() 

def quadrant_rate_plot(zenith, azimuth, n_files, graph, full_area):
    """
    Bin reconstructed tracks into NE/NW/SW/SE azimuth quadrants,
    compute rate (Hz), and plot as a bar chart.
    """
    if full_area:
        name_add_on = "Full Area"
    else:
        name_add_on = "Fiducial"

    # Define quadrants by azimuth range (degrees, 0=East, CCW)
    # Your track_angles uses arctan2(dy, dx) % 360, so 0=East, 90=North

    azimuth = (azimuth) % 360
    sec_in_day = 24 * 3600
    quadrant_masks = {
        'NE': (azimuth >= 0)   & (azimuth < 90),
        'NW': (azimuth >= 90)  & (azimuth < 180),
        'SW': (azimuth >= 180) & (azimuth < 270),
        'SE': (azimuth >= 270) & (azimuth < 360),
    }

    total_seconds = n_files * RUN_SECONDS
    quad_names  = list(quadrant_masks.keys())
    quad_counts = [np.sum(quadrant_masks[q]) for q in quad_names]
    quad_rates  = [(c / total_seconds) for c in quad_counts]
    quad_errors = [(np.sqrt(c) / total_seconds) for c in quad_counts]  # Poisson sqrt(N)

    fig, ax = plt.subplots(figsize=(8, 6))
    x_pos = np.arange(len(quad_names))

    ax.bar(x_pos, quad_rates, yerr=quad_errors, capsize=6,
           color='steelblue', alpha=0.8, label='Observed rate')

    ax.set_xticks(x_pos)
    ax.set_xticklabels(quad_names, fontsize=14)
    ax.set_xlabel('Quadrant', fontsize=14)
    ax.set_ylabel('Rate (Hz)', fontsize=14)
    ax.set_title(f'Quadrant Muon Rate {name_add_on}', fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    # Right axis: total counts
    ax2 = ax.twinx()
    ax2.set_ylim([y * total_seconds for y in ax.get_ylim()])
    ax2.set_ylabel('Total counts', fontsize=14)

    filename = f"{RUN_NAME}_Quadrant_Rate_{name_add_on}"
    filepath = os.path.join(SAVE_FOLDER, filename)
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    print(f'{filename} saved')

    if graph:
        plt.show()
    plt.clf()
    plt.close(fig)

def layer_hit_heatmap(data, graph, full_area):

    if full_area:
        BAR_EDGES = np.arange(8.25, 1064.25+16.5, 16.5)     
        h_line_min = 16.5
        h_line_max = 1056.0
        v_line_min = 16.5
        v_line_max = 1056.0
        name_add_on = "Full Area"
    else:
        BAR_EDGES = np.arange(8.25+16.5, 1064.25, 16.5)     
        h_line_min = 33
        h_line_max = 1039.5
        v_line_min = 33
        v_line_max = 1039.5
        name_add_on = "Fiducial"

    tracks_top = np.array([t[0] for t in data])
    x_top = tracks_top[:, 0]
    y_top = tracks_top[:, 1]
    #print(f"x range: {x_top.min()} to {x_top.max()}")  # expecting 0, 1023 or similar
    #print(f"Unique x values: {len(np.unique(x_top))}")
    #print(f"y range: {y_top.min()} to {y_top.max()}")  # expecting 0, 1023 or similar
    #print(f"Unique y values: {len(np.unique(y_top))}")
    
    plt.hist2d(x_top, y_top, bins=BAR_EDGES, cmap="viridis")
    
    # plt.axhline(h_line_min, color='red', linestyle='--', label='Horizontal Line')
    # plt.axhline(h_line_max, color='red', linestyle='--', label='Horizontal Line')
    # plt.axvline(v_line_min, color='red', linestyle='--', label='Vertical Line')
    # plt.axvline(v_line_max, color='red', linestyle='--', label='Vertical Line')

    plt.colorbar(label="count")
    plt.title(f'Top Layer Hit Positions {name_add_on}')

    filename = f"{RUN_NAME}_Top_xy_Heatmap_{name_add_on}"
    print(f'{filename} saved')
    filepath = os.path.join(SAVE_FOLDER, filename)
    plt.savefig(filepath, dpi=300)

    if graph:
        plt.show()
    
    plt.clf() 

    tracks_bottom = np.array([t[1] for t in data])
    x_bottom = tracks_bottom[:, 0]
    y_bottom = tracks_bottom[:, 1]
    print(f"x range: {x_bottom.min()} to {x_bottom.max()}")  # expecting 0, 1023 or similar
    print(f"Unique x values: {len(np.unique(x_bottom))}")
    print(f"y range: {y_bottom.min()} to {y_bottom.max()}")  # expecting 0, 1023 or similar
    print(f"Unique y values: {len(np.unique(y_bottom))}")

    plt.hist2d(x_bottom, y_bottom, bins=BAR_EDGES, cmap="viridis")
    # plt.axhline(h_line_min, color='red', linestyle='--', label='Horizontal Line')
    # plt.axhline(h_line_max, color='red', linestyle='--', label='Horizontal Line')
    # plt.axvline(v_line_min, color='red', linestyle='--', label='Vertical Line')
    # plt.axvline(v_line_max, color='red', linestyle='--', label='Vertical Line')
    plt.colorbar(label="count")
    plt.title(f'Bottom Layer Hit Positions {name_add_on}')

    filename = f"{RUN_NAME}_Bottom_xy_Heatmap_{name_add_on}"
    print(f'{filename} saved')
    filepath = os.path.join(SAVE_FOLDER, filename)
    plt.savefig(filepath, dpi=300)

    if graph:
        plt.show()
    
    plt.clf() 

def main(config):

    coincidence_files = folder_reader(SUB_DATA_FOLDER_PATH, file_max = 10, file_size_in_MB=None)
    print(f"Found {len(coincidence_files)} files for run {RUN_NAME}.")
    events    = read_coincidence_file(coincidence_files)
    ch_to_bar = build_ch_to_bar(config.LAYER_MAPS, config.OFFSETS)
    decoded   = decode_events(events, ch_to_bar)

    print(f"Total PETsys coincidence events: {len(events)}")
    print(f"Raw PETsys Coincidence rate: {len(events) / (len(coincidence_files) * RUN_SECONDS)} Hz")

    cuts = FilteredData(decoded,
                        max_bars_per_layer=5,
                        max_time_ns=20,
                        plot_time_cut=True,
                        save_folder=SAVE_FOLDER,
                        run_name=RUN_NAME,
                        show_plots=False)
    print(cuts.summary())

    bar_hits_per_layer(cuts.blocks, True, graph = False)
    
    bar_frequency_and_qdc_distribution_plots(cuts.blocks, graph= False)

    tracks_full = build_tracks(cuts.blocks, accept=bars_adjacent, cut_edge_bars=False)
    tracks_fiducial = fiducial_filter(tracks_full)

    zenith_full, azimuth_full = compute_angle_distributions(tracks_full)
    zenith_fiducial, azimuth_fiducial = compute_angle_distributions(tracks_fiducial)
    
    zenith_cut = 5
    vertical_muons_cut = zenith_fiducial[zenith_fiducial <= zenith_cut]
    print(f'Zenith Cut: {(len(vertical_muons_cut))}')
    vertical_muon_flux = len(vertical_muons_cut) / ((len(coincidence_files) * RUN_SECONDS) * (np.pi * (FIDUCIAL_RADIUS/1000)**2) * (2*np.pi * (1 - np.cos(np.radians(zenith_cut)))))
    print(f'Vertical Muon Flux: {vertical_muon_flux:.6e} muons/m^2/s/sr')
    
    azimuth_histogram = np.histogram(azimuth_fiducial, bins=45)
    print(f"Azimuthal Mean Variation: {np.std(azimuth_histogram[0])/np.mean(azimuth_histogram[0])*100:.2f}%")
    print(f"Maximum Azimuthal Variation: {(np.max(azimuth_histogram[0])-np.min(azimuth_histogram[0]))/np.min(azimuth_histogram[0])*100:.2f}%")

    print(f"Full Area Coincidence Rate: {len(tracks_full)/ (len(coincidence_files) * RUN_SECONDS)} Hz")
    print(f"Fiducial Area Coincidence Rate: {len(tracks_fiducial)/ (len(coincidence_files) * RUN_SECONDS)} Hz")

    #anglular_heatmap(zenith_fiducial, azimuth_fiducial, graph=False, full_area=False)
    #anglular_heatmap(zenith_full, azimuth_full, graph=False, full_area=True)

    #quadrant_rate_plot(zenith_fiducial, azimuth_fiducial, n_files=len(coincidence_files), graph=False, full_area=False)
    #quadrant_rate_plot(zenith_full, azimuth_full, n_files=len(coincidence_files), graph=False, full_area=True)

    zenith_and_azimuth_distribution_plot(zenith_fiducial, azimuth_fiducial, graph= False, full_area=False)
    zenith_and_azimuth_distribution_plot(zenith_full, azimuth_full, graph= False, full_area=True)

    layer_hit_heatmap(tracks_fiducial, graph=False, full_area=False)
    layer_hit_heatmap(tracks_full, graph=False, full_area=True)

    write_to_stats_file(cuts, len(events), len(coincidence_files), tracks_full, tracks_fiducial)

if __name__ == '__main__':
    #FILE_PATH   = 'KNVA-20260514-01-00079_coinc.dat'

    DETECTOR    = 'KN2'   # 'KN1' or 'KN2' -> loads <DETECTOR>_initialization.py
    config      = load_detector_config(DETECTOR)

    DATA_FOLDER_PATH = r"C:\\Users\\AlexClark\\Documents\\KoNova\\PETsys_Data"
    SAVE_FOLDER_PATH = r"C:\\Users\\AlexClark\\Documents\\KoNova\\PETsys_Plots"
    RUN_NAME    = 'KN2_Blue_Sky'
    SAVE_RUN_NAME = 'KN2_Blue_Sky_090526' #RUN_NAME
    SUB_DATA_FOLDER_PATH = os.path.join(DATA_FOLDER_PATH, RUN_NAME)
    SAVE_FOLDER = os.path.join(SAVE_FOLDER_PATH, SAVE_RUN_NAME)
    os.makedirs(SAVE_FOLDER, exist_ok=True) 
    print(f'Save directory created {SAVE_FOLDER}')

    main(config)
