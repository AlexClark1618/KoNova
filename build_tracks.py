import numpy as np
import time
import os
import matplotlib.pyplot as plt

# ----------------------------- Initial Variables -----------------------------
BAR_CH_MAP  = 'bar_ch_map.csv'
OFFSETS     = {1: 64, 2: 320, 3: 576, 4: 832} #Layer-ch num offset dict, e.g. layer 1 ch 0-63, layer 2 ch 320-383, etc.
BAR_SEP     = 16.5                       # mm, bar centre-to-centre
DELTA_Z     = 11.75 * 2.54 * 10          # mm, x-top to x-bottom (measured 4/29/2026)
RUN_SECONDS = 600
N_BARS      = 64
delta_x = 0.01 *10 #Dan's measured shifts between layers and stretch in bars
delta_y = 0.35 *10 
delta_x_stretch = 1.003 
delta_y_stretch = 0.9965 


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

def build_ch_to_bar(map_path, offsets):
    """Map scaled channel number -> (layer, bar)."""
    ch_to_bar = {}
    with open(map_path, 'r') as f:
        next(f)  # skip header
        for line in f:
            bar, base_ch = map(int, line.strip().split(','))
            for layer, off in offsets.items():
                ch_to_bar[base_ch + off] = (layer, bar)
    return ch_to_bar


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
def per_layer_stats(decoded):
    """Return (freq, qdc) dicts keyed by layer."""
    freq_dist = {l: {b: 0 for b in range(1, N_BARS + 1)} for l in range(1, 5)}
    qdc_dist  = {l: [] for l in range(1, 5)}
    #hits_by_layer = {l: [h for h in block if h[0] == l] for l in range(1, 5)}

    for block in decoded:
        for layer, bar, qdc_value, _ in block:
            freq_dist[layer][bar] += 1
            qdc_dist[layer].append(qdc_value)

    return freq_dist, qdc_dist

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

    if len(set(bars)) > 1: 
        #If multiple bars hit, return random position between lowest and highest hit bars hit
        low  = float(bars.min()* sep)   # Center of lowest hit bar
        high = float(bars.max() * sep)  # Center of highest hit bar
        bar_pos = np.random.uniform(low, high)
        return bar_pos
    
    else:
        return int(bars) * sep # Bar 1 center = 16.5, Between Bar 1 and 2 = 24.75, ... Bar 64 center = 1056.0


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

'''
def hit_count_per_layer(decoded):
    for block in decoded:
        print(decoded)
        hits_by_layer = {l: [h for h in block if h[0] == l] for l in range(1, 5)}
'''

def build_tracks(decoded, accept=bars_adjacent, cut_edge_bars=True):
    """
    Build [(x_top, y_top), (x_bottom, y_bottom)] per event.

    accept : optional callable(layer_hits) -> bool, applied per layer.
             Event is dropped if any required layer fails or is empty.
    """
    tracks = []
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
        tracks.append([top, bottom])
    print(f"Events skipped (empty layers): {skipped_none}")
    print(f"Events skipped (non-adjacent bars): {skipped_adjacency}")
    print(f"Events skipped (edge bars): {skipped_edge}")
    print(f"Events kept: {len(tracks)}")
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
    angles = np.array([track_angles(t, b, delta_z) for t, b in tracks])
    return angles[:, 0], angles[:, 1]   # zenith, azimuth


def background_subtraction(data1, data2):
    BAR_EDGES = np.arange(8.25, 1064.25+16.5, 16.5)     
    edges = np.arange(0,1072+16.5, 16.5)
    print(len(BAR_EDGES))
    #Signal
    tracks_top = np.array([t[0] for t in data1])
    x_top = tracks_top[:, 0]
    y_top = tracks_top[:, 1]

    signal_top_hist, x_edges, y_edges = np.histogram2d(x_top, y_top, bins = BAR_EDGES)
    print(x_edges)
    print(y_edges)
    tracks_bottom = np.array([t[1] for t in data1])
    x_bottom = tracks_bottom[:, 0]
    y_bottom = tracks_bottom[:, 1]

    signal_bottom_hist, _ , _ = np.histogram2d(x_bottom, y_bottom, bins = BAR_EDGES)

    #Background
    tracks_top = np.array([t[0] for t in data2])
    x_top = tracks_top[:, 0]
    y_top = tracks_top[:, 1]

    background_top_hist, _, _ = np.histogram2d(x_top, y_top, bins = BAR_EDGES)

    tracks_bottom = np.array([t[1] for t in data2])
    x_bottom = tracks_bottom[:, 0]
    y_bottom = tracks_bottom[:, 1]
    
    background_bottom_hist, _, _ = np.histogram2d(x_bottom, y_bottom, bins = BAR_EDGES)

    #Subtraction
    unshielded = np.ones_like(signal_top_hist, dtype=bool)
    unshielded[400:600, 400:600] = False

    k = signal_top_hist[unshielded].sum() / background_top_hist[unshielded].sum()
    print(f"k = {k:.6f}")
    print(f"Net residual: {(signal_top_hist - k * background_top_hist).sum():.1f}")    
    top_residual = signal_top_hist - k* background_top_hist
    print("Top residual sum",np.sum(top_residual))


    var   = signal_top_hist + k**2 * background_top_hist                     # variance
    sigma = np.sqrt(var)                                   # sqrt element-wise
    Z_top = np.where(sigma > 0, top_residual / sigma, 0.0)

    fig, ax = plt.subplots()
    im = ax.pcolormesh(x_edges, y_edges, top_residual.T, cmap="RdBu_r")
    plt.colorbar(im, ax=ax, label="Flux difference (counts)")
    plt.xlim(0, 1072)
    plt.xlabel("x (mm)")
    plt.ylim(0, 1072)
    plt.xlabel("y (mm)")

    detector_length = 1072
    plt.axhline(y= detector_length-469.9, color='red', linestyle='--', label='Horizontal Line')
    plt.axhline(y = detector_length-673.1, color='red', linestyle='--', label='Horizontal Line')
    plt.axvline(x = detector_length-469.9, color='red', linestyle='--', label='Vertical Line')
    plt.axvline(x = detector_length-673.1, color='red', linestyle='--', label='Vertical Line')
    plt.show()

    fig, ax = plt.subplots()
    im = ax.pcolormesh(x_edges, y_edges, Z_top.T, cmap="RdBu_r")
    plt.colorbar(im, ax=ax, label="Z Map")
    plt.show()

    k = signal_bottom_hist.sum() / background_bottom_hist.sum()
    print(f"k = {k:.6f}")
    print(f"Net residual: {(signal_bottom_hist - k * background_bottom_hist).sum():.1f}")    
    bottom_residual = signal_bottom_hist - k* background_bottom_hist
    print("Bottom residual sum",np.sum(bottom_residual))

    var   = signal_bottom_hist + k**2 * background_bottom_hist                     # variance
    sigma = np.sqrt(var)                                   # sqrt element-wise
    Z_bottom = np.where(sigma > 0, bottom_residual / sigma, 0.0)

    fig, ax = plt.subplots()
    im = ax.pcolormesh(x_edges, y_edges, bottom_residual.T, cmap="RdBu_r")
    plt.colorbar(im, ax=ax, label="Flux difference (counts)")
    plt.show()

    fig, ax = plt.subplots()
    im = ax.pcolormesh(x_edges, y_edges, Z_bottom.T, cmap="RdBu_r")
    plt.colorbar(im, ax=ax, label="Z Map")
    plt.show()
# ----------------------------- Plots -----------------------------

def layer_hit_heatmap(data, graph, full_area):

    if full_area:
        BAR_EDGES = np.arange(0, 1072.5+16.5, 16.5)
        h_line_min = 16.5
        h_line_max = 1056.0
        v_line_min = 16.5
        v_line_max = 1056.0
        name_add_on = "Full Area"
    else:
        BAR_EDGES = np.arange(16.5, 1056+16.5, 16.5)   # bar 1 right edge to bar 63 right edge
        h_line_min = 33
        h_line_max = 1039.5
        v_line_min = 33
        v_line_max = 1039.5
        name_add_on = "Fiducial"

    tracks_top = np.array([t[0] for t in data])
    x_top = tracks_top[:, 0]
    y_top = tracks_top[:, 1]
    print(f"x range: {x_top.min()} to {x_top.max()}")  # expecting 0, 1023 or similar
    print(f"Unique x values: {len(np.unique(x_top))}")
    print(f"y range: {y_top.min()} to {y_top.max()}")  # expecting 0, 1023 or similar
    print(f"Unique y values: {len(np.unique(y_top))}")
    
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

def main():

    SUB_DATA_FOLDER_PATH = os.path.join(DATA_FOLDER_PATH, RUN_NAME)
    coincidence_files = folder_reader(SUB_DATA_FOLDER_PATH, file_max=5)
    print(f"Found {len(coincidence_files)} files for run {RUN_NAME}.")
    events    = read_coincidence_file(coincidence_files)
    ch_to_bar = build_ch_to_bar(BAR_CH_MAP, OFFSETS)
    decoded   = decode_events(events, ch_to_bar)

    print(f"Total coincidence events: {len(events)}")
    print(f"Coincidence rate: {len(events) / (len(coincidence_files) * RUN_SECONDS)} Hz")

    tracks_full_lead = build_tracks(decoded, accept=bars_adjacent, cut_edge_bars=False)

    SUB_DATA_FOLDER_PATH = os.path.join(DATA_FOLDER_PATH, BACKGROUND)
    coincidence_files = folder_reader(SUB_DATA_FOLDER_PATH, file_max=5)
    print(f"Found {len(coincidence_files)} files for run {BACKGROUND}.")
    events    = read_coincidence_file(coincidence_files)
    ch_to_bar = build_ch_to_bar(BAR_CH_MAP, OFFSETS)
    decoded   = decode_events(events, ch_to_bar)

    print(f"Total coincidence events: {len(events)}")
    print(f"Coincidence rate: {len(events) / (len(coincidence_files) * RUN_SECONDS)} Hz")

    tracks_full_background = build_tracks(decoded, accept=bars_adjacent, cut_edge_bars=False)

    background_subtraction(tracks_full_lead, tracks_full_background)
    #tracks_fiducial = build_tracks(decoded, accept=bars_adjacent, cut_edge_bars=True)
    #tracks_fiducial = fiducial_filter(tracks_full)


    print(f"Full Area Coincidence Rate: {len(tracks_full)/ (len(coincidence_files) * RUN_SECONDS)} Hz")
    #print(f"Fiducial Area Coincidence Rate: {len(tracks_fiducial)/ (len(coincidence_files) * RUN_SECONDS)} Hz")

    #layer_hit_heatmap(tracks_fiducial, graph=False, full_area=False)
    layer_hit_heatmap(tracks_full, graph=False, full_area=True)

if __name__ == '__main__':
    FILE_PATH   = 'KNVA-20260514-01-00079_coinc.dat'
    DATA_FOLDER_PATH = r"C:\\Users\\aclark2\\Desktop\\KoNova Code\\PETsys Data"
    SAVE_FOLDER_PATH = r"C:\\Users\\aclark2\\Desktop\\KoNova Code\\PETsys Plots"
    RUN_NAME    = 'Lead'
    BACKGROUND = 'Loading Dock'
    SAVE_RUN_NAME = 'Practice' #RUN_NAME
    #SUB_DATA_FOLDER_PATH = os.path.join(DATA_FOLDER_PATH, RUN_NAME)
    SAVE_FOLDER = os.path.join(SAVE_FOLDER_PATH, RUN_NAME)
    os.makedirs(SAVE_FOLDER, exist_ok=True) 
    print(f'Save directory created {SAVE_FOLDER}')

    main()