# KoNova Muon Detection Analysis

A Python analysis pipeline for processing and visualizing cosmic ray muon data collected with a four layer Hodoscope detector (KoNova) using the PETsys DAQ system, fermilab triangular scintillator bars, kurary y-11 wavelength shifting fiber.

---

## Overview

The KoNova detector is a two-layer x-y scintillator bar system that reconstructs cosmic ray muon tracks in 3D. Each layer consists of 64 bars with a 16.5 mm center-to-center spacing. The detector records coincidence events — hits that occur across all four layers (x-top, y-top, x-bottom, y-bottom) within a time window — which are then decoded, filtered, and analyzed to reconstruct muon trajectories and spatial flux maps.

The primary application demonstrated here is **muon attenuation imaging**: comparing a signal run (with a lead castle shielding object placed on the detector) against a background run to produce statistically significant flux difference and z-score maps.

---

## Detector Layout

```
Top plane:    Layer 3 (x)  ×  Layer 4 (y)
                        ↕  DELTA_Z = 298.45 mm
Bottom plane: Layer 1 (x)  ×  Layer 2 (y)
```

- 64 bars per layer, 16.5 mm pitch → ~1056 mm active area
- Spatial resolution: **16.5 mm** (one bar width)
- Layer separation measured: 11.75 inches (4/29/2026)

---

## Repository Structure

```
├── KoNova_Analysis.py       # Main analysis script — angles, rates, heatmaps
├── build_tracks.py          # Background subtraction and z-score mapping
├── read_coin_file.py        # Early prototype reader and angle calculator
├── bar_ch_map.csv           # Bar number ↔ base channel mapping
└── README.md
```

---

## Pipeline

### 1. File Reading
`folder_reader()` + `read_coincidence_file()`

Reads all `_coinc.dat` files from a run directory. Each file contains coincidence event blocks — groups of hits (timestamp, QDC, channel) that fired within the DAQ coincidence window.

### 2. Channel Decoding
`build_ch_to_bar()` + `decode_events()`

Maps raw DAQ channel numbers to `(layer, bar)` pairs using `bar_ch_map.csv` and layer offsets:

| Layer | Channel Offset |
|-------|---------------|
| 1     | 64            |
| 2     | 320           |
| 3     | 576           |
| 4     | 832           |

### 3. Track Building
`build_tracks()`

Reconstructs `[(x_top, y_top), (x_bottom, y_bottom)]` track pairs per event. Events are filtered by:
- **Empty layer cut** — all four layers must have at least one hit
- **Adjacency cut** — hits within a layer must be on consecutive bars (no split clusters)
- **Edge bar cut** (optional) — drops events hitting bar 1 or bar 64

Hit position within a layer is taken as the bar center (`bar × 16.5 mm`). If multiple adjacent bars fire, position is drawn uniformly between the outermost bar centers.

### 4. Angle Reconstruction
`track_angles()` + `compute_angle_distributions()`

Computes **zenith** (from vertical) and **azimuth** (CCW from +x) angles from the top-to-bottom displacement vector, with small geometric corrections for inter-layer alignment shifts measured by Dan:

```python
delta_x = 0.1 mm,  delta_y = 3.5 mm
stretch_x = 1.003,  stretch_y = 0.9965
```

### 5. Background Subtraction & Z-Score Mapping
`background_subtraction()`

Compares signal (lead) and background runs pixel-by-pixel:

1. **Histogram** both runs into 2D spatial grids using `np.histogram2d` with 16.5 mm bins
2. **Compute normalization** `k = signal[unshielded].sum() / background[unshielded].sum()` using only the unshielded region outside the lead castle
3. **Compute residual** `residual = signal - k × background`
4. **Compute z-score** per pixel:

```
Z = (signal - k × background) / sqrt(signal + k² × background)
```

The denominator is the Poisson uncertainty on the residual, propagating uncertainties from both the signal and scaled background.

> **Note on masking:** The shielded region (lead castle footprint: 431.8–635 mm) is defined in **bin indices**, not mm coordinates. Converting: `bin = int(mm / 16.5)` gives bins 26–38.

### 6. Outputs

| Plot | Description |
|------|-------------|
| Flux Difference Map | Raw `signal - k×background` counts per pixel |
| Z-Score Map | Per-pixel statistical significance |
| Z-Score Distribution (Shielded) | Should be Gaussian shifted negative if lead is detected |
| Z-Score Distribution (Unshielded) | Should be Gaussian centered on 0 if normalization is correct |
| Bar Hit Frequency | Per-layer bar hit rates for detector diagnostics |
| QDC Distribution | Charge distribution per layer |
| Zenith/Azimuth Distributions | Angular flux distributions |
| Angular Heatmap (Polar) | Raw and solid-angle-corrected angular flux |
| Quadrant Rate Plot | N/S/E/W muon rate asymmetry |
| Layer Hit Heatmaps | 2D spatial hit maps for signal and background |

---

## Configuration

Set these variables in the `if __name__ == '__main__':` block before running:

```python
DATA_FOLDER_PATH = r"path\to\PETsys Data"   # Root data directory
SAVE_FOLDER_PATH = r"path\to\PETsys Plots"  # Root output directory
RUN_NAME         = 'Lead'                    # Signal run subfolder
BACKGROUND       = 'Loading Dock'            # Background run subfolder
SAVE_RUN_NAME    = 'Practice'               # Output subfolder name
RUN_SECONDS      = 600                       # Duration of each DAQ file (s)
```

---

## Requirements

```
numpy
matplotlib
```

Standard library: `os`, `collections`, `time`

---

## Key Physical Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `BAR_SEP` | 16.5 mm | Bar center-to-center pitch |
| `DELTA_Z` | 298.45 mm | Vertical separation between top and bottom planes |
| `N_BARS` | 64 | Bars per layer |
| `OFFSETS` | {1:64, 2:320, 3:576, 4:832} | DAQ channel offsets per layer |
| Spatial resolution | 16.5 mm | Limited by bar width |
| Lead castle region | 431.8–635 mm | Both x and y (bins 26–38) |
