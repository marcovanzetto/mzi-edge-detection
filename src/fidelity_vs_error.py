"""
Fidelity of the uncompensated circuit as a function of the fabrication error on n_eff.

The parameters, the input image and the functions are the ones of mzi_edge_detection.py:
only the error on n_eff is changed, while the errors on DL and f_c are kept.

Usage: python src/fidelity_vs_error.py [path/to/image]
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq, fftshift

import mzi_edge_detection as mzi


# =============================================================================
# SETTINGS
# =============================================================================

is_plot = True # to show the plot of F as a function of the error on n_eff
dn_plot_max = 0.015 # the plot covers the errors in [-dn_plot_max, dn_plot_max]
dn_ranges = [0.0005, 0.001, 0.003, 0.006, 0.012, 0.024] # half widths of the uniform distributions of the error
n_samples = 20000 # number of random errors for each distribution
seed = 0 # for reproducible results

output_path = Path(__file__).resolve().parent.parent / "results" / "fidelity_vs_error.png"


# =============================================================================
# SIGNAL BANDWIDTH
# =============================================================================

# same procedure of compute_derivative: maximum bandwidth over the rows of the image

dim_col = mzi.img_array.shape[1]
dt = (dim_col / mzi.f_s) / (dim_col - 1)
f = fftfreq(dim_col, d=dt) + mzi.f_c

BW_max = max(
    mzi.bandwidth_threshold(fftshift(f), fftshift(np.abs(fft(row / np.sqrt(2)))))
    for row in mzi.img_array
)


# =============================================================================
# FIDELITY
# =============================================================================

f_rel_dense = np.linspace(-BW_max/2, BW_max/2, int(1e4))
f_dense = f_rel_dense + mzi.f_c

H_id_dense = 1j * mzi.DphiL_id(f_rel_dense) / 2
norm_H_id = mzi.trapezoid(np.abs(H_id_dense)**2, x=f_dense)

# same overlap integral of compute_derivative, for an array of errors on n_eff and without compensation (dT = 0)
def fidelity_uncompensated(dn_eff):

    dn_eff = np.atleast_1d(dn_eff)[:, None]

    DphiL = (2 * np.pi * f_dense / mzi.c) * (mzi.n_eff_id + dn_eff) * mzi.DL
    H_MZI_dense = 1j * np.exp(-1j * DphiL / 2) * np.sin(DphiL / 2)

    overlap = mzi.trapezoid(np.conj(H_id_dense) * H_MZI_dense, x=f_dense, axis=1)
    norm_H = mzi.trapezoid(np.abs(H_MZI_dense)**2, x=f_dense, axis=1)

    return np.abs(overlap)**2 / (norm_H_id * norm_H)

# the errors are processed in blocks to limit the memory usage
def fidelity_in_blocks(dn_eff, block=200):
    n_blocks = max(1, len(dn_eff) // block)
    return np.concatenate([fidelity_uncompensated(b) for b in np.array_split(dn_eff, n_blocks)])


# =============================================================================
# RANDOM ERRORS
# =============================================================================

rng = np.random.default_rng(seed)

print("===== Uncompensated circuit with random errors on n_eff ========")
print(f"BW_max = {BW_max:.3g} Hz")
print(f"period of F in dn_eff = n_eff_id/m = {mzi.n_eff_id/mzi.m:.4g}")
print(f"{n_samples} errors drawn uniformly for each range\n")
print("range of dn_eff    F median    F > 0.1    F > 0.9")

for dn_range in dn_ranges:
    F = fidelity_in_blocks(rng.uniform(-dn_range, dn_range, n_samples))
    print(f"+/- {dn_range:<10g} {np.median(F):11.2g} {100*np.mean(F > 0.1):9.1f}% {100*np.mean(F > 0.9):9.1f}%")


# =============================================================================
# PLOT
# =============================================================================

dn_axis = np.linspace(-dn_plot_max, dn_plot_max, 3001)
F_axis = fidelity_in_blocks(dn_axis)

plt.figure(figsize=(8, 4))
plt.plot(dn_axis, F_axis, color='b')
plt.plot(mzi.dn_eff_err, fidelity_uncompensated(mzi.dn_eff_err)[0], 'o', color='r',
         label=f'default value ($\\delta n_\\mathrm{{eff}}$ = {mzi.dn_eff_err:g})')
plt.xlabel('$\\delta n_\\mathrm{eff}$')
plt.ylabel('$F$ (uncompensated)')
plt.xlim(-dn_plot_max, dn_plot_max)
plt.ylim(0, 1.05)
plt.legend(fontsize=11)
plt.grid()
plt.tight_layout()
plt.savefig(output_path, dpi=110)

if is_plot:
    plt.show()
