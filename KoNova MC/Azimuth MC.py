"""
Monte Carlo: misalignment-induced azimuthal modulation in a 2-layer circular
hodoscope. Handles BOTH tilts and relative layer shifts.

Model
-----
Two circular fiducial layers of radius r_fid, separated in z by dz.
A track is defined by its true direction (phi_true, theta_true) and a vertex
position on the detector face. Its horizontal displacement between the two
layers is:

    rho   = dz * tan(theta)            (horizontal travel between layers)
    dx    = rho * cos(phi)
    dy    = rho * sin(phi)

Detector misalignment of layer 2 relative to layer 1 combines two effects:

  TILT (offset proportional to dz):
    alpha  = rotation about the x-axis (pitch)  -> offsets dy by -alpha*dz
    beta   = rotation about the y-axis (yaw)    -> offsets dx by +beta *dz

  SHIFT (fixed translation, independent of dz):
    delta_x = layer-2 offset in x   -> offsets dx by +delta_x
    delta_y = layer-2 offset in y   -> offsets dy by +delta_y

So the measured hit differences are:
    dx' = dx + beta*dz + delta_x
    dy' = dy - alpha*dz + delta_y

Reconstructed azimuth:
    phi_meas = atan2(dy', dx')

Key distinction: tilt offsets scale with dz, shift offsets do not. Both give a
cot(theta) zenith dependence, but only the shift scales as 1/dz -- that's the
lever to separate them.

Acceptance: both hits must fall inside their fiducial circles. The vertex (hit1)
is placed uniformly on the fiducial disk; hit2 = hit1 + true displacement. Both
must lie within r_fid.

Outputs a histogram of reconstructed azimuth (ideal vs misaligned) plus a
simultaneous dipole fit reporting the cos and sin components.
"""

import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Detector / simulation parameters
# ----------------------------------------------------------------------
R_FID   = 0.5          # fiducial radius (m)
DZ      = 0.30        # layer separation in z (m)
ALPHA   = np.deg2rad(0.0)   # pitch: rotation about x-axis (rad)
BETA    = np.deg2rad(1.0)   # yaw:   rotation about y-axis (rad)
DELTA_X = 0#-0.0018#-0.0015         # relative layer-2 shift in x (m)
DELTA_Y = 0#0.0028#0.0025         # relative layer-2 shift in y (m)
THETA_MAX = np.deg2rad(70.0)  # max zenith angle (rad)
N_TRACKS  = 2_000_000  # number of trial tracks
N_BINS    = 72         # azimuth histogram bins
RNG_SEED  = 12345

rng = np.random.default_rng(RNG_SEED)


def sample_zenith(n):
    """Rejection-sample zenith from the distribution  f(theta) ~ cos^3(theta)*sin(theta)
    on [0, pi/2), then keep only theta <= THETA_MAX.

    This is the cosmic-ray-like angular distribution (cos^2 flux times the
    cos*sin solid-angle/projection factor). We draw candidates, accept with
    probability proportional to cos^3*sin (peak value of that function is used
    to normalise the acceptance so it stays a valid probability <= 1)."""
    # peak of cos^3(theta)*sin(theta) occurs at theta = atan(1/sqrt(3)) = 30 deg
    theta_peak = np.arctan(1.0 / np.sqrt(3.0))
    f_max = np.cos(theta_peak)**3 * np.sin(theta_peak)

    out = np.empty(0, dtype=float)
    while out.size < n:
        cand = rng.uniform(0.0, np.pi / 2.0, size=n)
        prob = (np.cos(cand)**3) * np.sin(cand) / f_max
        keep = rng.uniform(0.0, 1.0, size=n) < prob
        cand = cand[keep]
        cand = cand[cand <= THETA_MAX]          # enforce detector max zenith
        out = np.concatenate([out, cand])
    return out[:n]


def sample_tracks(n):
    """Sample truth tracks: uniform phi, cos^3*sin zenith, uniform vertex."""
    # azimuth: uniform on [0, 2pi)
    phi = rng.uniform(0.0, 2.0 * np.pi, size=n)

    # zenith: cos^3(theta)*sin(theta) via rejection sampling
    theta = sample_zenith(n)

    # vertex position of hit 1: uniform over the fiducial disk of layer 1
    # (sample uniformly in the disk via sqrt trick)
    r1 = R_FID * np.sqrt(rng.uniform(0.0, 1.0, size=n))
    a1 = rng.uniform(0.0, 2.0 * np.pi, size=n)
    x1 = r1 * np.cos(a1)
    y1 = r1 * np.sin(a1)

    return phi, theta, x1, y1


def reconstruct(phi, theta, x1, y1, alpha, beta, delta_x, delta_y):
    """Apply tilt + shift, compute hit2, apply acceptance, return reconstructed phi."""
    rho = DZ * np.tan(theta)          # horizontal displacement magnitude
    dx = rho * np.cos(phi)            # true displacement
    dy = rho * np.sin(phi)

    # hit 2 absolute position (before misalignment) = hit1 + true displacement
    x2 = x1 + dx
    y2 = y1 + dy

    # ---- acceptance: both hits inside fiducial circles ----
    inside1 = (x1**2 + y1**2) <= R_FID**2
    inside2 = (x2**2 + y2**2) <= R_FID**2
    accept = inside1 & inside2

    # ---- misalignment on the measured differences ----
    #   tilt term  scales with dz:  +beta*dz (x),  -alpha*dz (y)
    #   shift term is fixed:        +delta_x (x),  +delta_y (y)
    dx_meas = dx + beta * DZ + delta_x
    dy_meas = dy - alpha * DZ + delta_y

    phi_meas = np.arctan2(dy_meas, dx_meas)
    phi_meas = np.mod(phi_meas, 2.0 * np.pi)   # wrap to [0, 2pi)

    return phi_meas, accept, theta


def main():
    phi, theta, x1, y1 = sample_tracks(N_TRACKS)

    # Ideal (no misalignment)
    phi_ideal, acc_ideal, _ = reconstruct(phi, theta, x1, y1,
                                          0.0, 0.0, 0.0, 0.0)
    # Misaligned (tilt + shift)
    phi_mis, acc_mis, _ = reconstruct(phi, theta, x1, y1,
                                      ALPHA, BETA, DELTA_X, DELTA_Y)

    phi_ideal_acc = phi_ideal[acc_ideal]
    phi_mis_acc   = phi_mis[acc_mis]

    bins = np.linspace(0.0, 2.0 * np.pi, N_BINS + 1)
    centers = 0.5 * (bins[:-1] + bins[1:])

    h_ideal, _ = np.histogram(phi_ideal_acc, bins=bins)
    h_mis,   _ = np.histogram(phi_mis_acc,   bins=bins)

    # peak-to-valley fractional modulation of the misaligned histogram
    ptv = (h_mis.max() - h_mis.min()) / (h_mis.max() + h_mis.min())

    # ---- plot ----
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))

    label_mis = (f"tilt a={np.rad2deg(ALPHA):.1f} b={np.rad2deg(BETA):.1f} deg,  "
                 f"shift dx={DELTA_X*100:.1f} dy={DELTA_Y*100:.1f} cm")

    ax[0].step(np.rad2deg(centers), h_ideal, where="mid",
               label="no misalignment", color="tab:gray")
    ax[0].step(np.rad2deg(centers), h_mis, where="mid",
               label=label_mis, color="tab:red")
    ax[0].set_xlabel("reconstructed azimuth (deg)")
    ax[0].set_ylabel("counts")
    ax[0].set_title(f"Azimuth distribution  (peak-to-valley = {ptv*100:.1f}%)")
    ax[0].legend(fontsize=8)
    ax[0].set_xlim(0, 360)

    # ratio to show the modulation shape cleanly
    ratio = h_mis / np.clip(h_ideal, 1, None)
    ax[1].step(np.rad2deg(centers), ratio, where="mid", color="tab:blue")
    ax[1].axhline(1.0, ls="--", color="k", lw=0.8)
    ax[1].set_xlabel("reconstructed azimuth (deg)")
    ax[1].set_ylabel("misaligned / ideal")
    ax[1].set_title("Modulation shape (dipole expected)")
    ax[1].set_xlim(0, 360)

    fig.tight_layout()
    fig.savefig("azimuth_tilt_mc.png", dpi=130)

    # ---- console summary ----
    print(f"tracks generated         : {N_TRACKS:,}")
    print(f"accepted (ideal)         : {acc_ideal.sum():,}")
    print(f"accepted (misaligned)    : {acc_mis.sum():,}")
    print(f"pitch alpha (deg)        : {np.rad2deg(ALPHA):.2f}")
    print(f"yaw   beta  (deg)        : {np.rad2deg(BETA):.2f}")
    print(f"shift delta_x (cm)       : {DELTA_X*100:.2f}")
    print(f"shift delta_y (cm)       : {DELTA_Y*100:.2f}")
    print(f"peak-to-valley mod       : {ptv*100:.1f}%")

    # quick dipole fit: N(phi) = A + B*cos(phi) + C*sin(phi)
    M = np.column_stack([np.ones_like(centers),
                         np.cos(centers), np.sin(centers)])
    coef, *_ = np.linalg.lstsq(M, h_mis, rcond=None)
    A, B, C = coef
    dipole_amp = np.hypot(B, C) / A
    dipole_phase = np.rad2deg(np.arctan2(C, B))
    print(f"fitted dipole cos coeff  : {B/A*100:+.1f}% of mean  (x-shift / pitch)")
    print(f"fitted dipole sin coeff  : {C/A*100:+.1f}% of mean  (y-shift / yaw)")
    print(f"fitted dipole amplitude  : {dipole_amp*100:.1f}% of mean")
    print(f"fitted dipole phase      : {dipole_phase:.1f} deg  "
          f"(peak of dN/dphi)")

    # display the interactive window (blocks until you close it)
    plt.show()


if __name__ == "__main__":
    main()