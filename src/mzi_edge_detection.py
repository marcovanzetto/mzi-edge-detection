import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft, fftfreq, fftshift
from PIL import Image

# np.trapz was renamed np.trapezoid in numpy 2.0, this works with both versions
trapezoid = np.trapezoid if hasattr(np, "trapezoid") else np.trapz


# =============================================================================
# TOGGLES
# =============================================================================

is_image_plot = True # to show the final image with highlighted edges
is_control_plots = False # to show the control plots for a row (or column, see last lines of this code)
is_compensated = True # to compensate the circuit with thermal phase shifters
is_ideal = False # to remove all the non idealities

nonideal = 0 if is_ideal else 1 # factor that multiplies the tolerances (0 in the ideal case)


# =============================================================================
# PHYSICAL PARAMETERS
# =============================================================================

c = 3e8 # vacuum light velocity

K_th = 1.86e-4 # thermo-optic coefficient

# refractive index
n_eff_id = 2.2
dn_eff_dw = 3e-3 # variation on n_eff due to width change
dn_eff_dh = 4.2e-3 # variation on n_eff due to height change
# fabrication error on n_eff due to tolerances in thickness and width of waveguides
# without compensation this value leaves the nearest zero about 9 GHz from the carrier, see README
dn_eff_err = 0.01194
dn_eff_wh = dn_eff_err*nonideal
n_eff = n_eff_id + dn_eff_wh

# carrier frequency
f_c_id = 1.93e14
df_c = 1e8*nonideal # tolerance
f_c = f_c_id + df_c

f_s = 1e11 # sampling frequency


# =============================================================================
# MZI
# =============================================================================

dDL = 1e-8*nonideal # tolerance of the path unbalance between the two arms of the MZI

dm_id = 1 # target of dm
m = 180 # order of the minimum of the TF

if dn_eff_wh <= 0: # the target minimum changes side with negative tolerance (with no error the m-th minimum is already aligned)
    dm_id = 0

DL_id = m * c / (n_eff_id * f_c_id) # we choose DL through m to align the real notch to the ideal one
DL = DL_id + dDL

FSR = c / (n_eff * DL)

Dm_n_uncomp = m*dn_eff_wh/n_eff_id # misalignment of the notches due only to dn (expressed in dm)
df = np.abs(f_c - m*FSR) # misalignment of the notches (expressed in frequency)
Dm_uncomp = df/FSR # misalignment of the notches (expressed in dm)

# ----------------
# COMPENSATION (only for tolerances on n_eff and DL, not on f_c)
# ----------------

# calculation of the dT needed to compensate at the target dm
dT = (dm_id * c * n_eff_id - m * c * dn_eff_wh - dDL * f_c_id * (n_eff_id + dn_eff_wh) * n_eff_id) / (m * c * K_th + n_eff_id * dDL * f_c_id * K_th)

if not is_compensated:
    dT = 0

dn_eff_T = K_th * dT
dn_eff = dn_eff_wh + dn_eff_T
n_eff += dn_eff_T

Dm_comp = (m*c*dn_eff+dDL*f_c*n_eff*n_eff_id)/(c*n_eff_id) # compensated dm (should be dm_id)

# for the TF
DphiL_id = lambda f: (2 * np.pi * f / c) * n_eff_id * DL_id
DphiL = lambda f: (2 * np.pi * f / c) * n_eff * DL

FSR = c / (n_eff * DL)

# detuning of the zero closest to the carrier (the zeros of the TF are at the integer multiples of FSR):
# after compensation it is the target zero (m + dm_id), without compensation it can be the (m-1)-th,
# the m-th or the (m+1)-th zero, depending on the sign and size of the error
df = np.abs(f_c - np.round(f_c / FSR) * FSR)


# =============================================================================
# IMAGE LOADING
# =============================================================================

# the image can be passed as command line argument, otherwise the default one in the images folder is used
default_image_path = Path(__file__).resolve().parent.parent / "images" / "iss.jpg"
image_path = sys.argv[1] if len(sys.argv) > 1 else default_image_path

image = Image.open(image_path)

# conversion of gray scale or palette images to RGB, so that the red color can be applied to the edges
if image.mode not in ("RGB", "RGBA"):
    image = image.convert("RGB")

original_image = np.array(image, dtype=float) / 255.0

img_array = np.array(
    image.convert("L"),
    dtype=float
) / 255.0  # conversion in gray tones

# adaptation of png format or jpg (transparency channel needed in the former case)
if image.mode == 'RGBA':
    red_color = np.array([1.0, 0.0, 0.0, 1.0])  # RGBA for .png
else:
    red_color = np.array([1.0, 0.0, 0.0]) # for .jpg


# =============================================================================
# BANDWIDTH ESTIMATION
# =============================================================================

# width of the band around the peak of the spectrum in which the amplitude stays above fraction times the peak
# (the band ends at the first frequency, on each side, at which the amplitude falls below this threshold)
def bandwidth_threshold(f, Y, fraction=1/100):

    Y = np.asarray(Y)
    f = np.asarray(f)

    # only positive frequencies (f is usually the absolute frequency f_rel + f_c, so all the points are kept)
    Y = Y[f>0]
    f = f[f>0]

    i0 = np.argmax(Y)

    threshold = fraction * Y[i0]

    # ----------------
    # left side
    # ----------------

    i = i0
    while i > 0 and Y[i] >= threshold:
        i -= 1

    if i == 0:
        f_left = f[0]
    else:
        # linear interpolation
        alpha = (threshold - Y[i]) / (Y[i+1] - Y[i])
        f_left = f[i] + alpha * (f[i+1]-f[i])

    # ----------------
    # right side
    # ----------------

    i = i0
    while i < len(Y)-1 and Y[i] >= threshold:
        i += 1

    if i == len(Y)-1:
        f_right = f[-1]
    else:
        alpha = (threshold - Y[i]) / (Y[i-1]-Y[i])
        f_right = f[i] + alpha * (f[i-1]-f[i])

    BW = f_right - f_left

    return BW


# =============================================================================
# DERIVATIVE CALCULATION
# =============================================================================

def compute_derivative(img_array, are_rows=True): # are_rows is False for the columns (transposed image)

    dim_row = img_array.shape[0]
    dim_col = img_array.shape[1]

    # ----------------
    # Temporal domain
    # ----------------

    T = dim_col / f_s # total time of the waveform to differentiate

    dt = T / (dim_col - 1)
    t = np.arange(dim_col) * dt

    # ----------------
    # Frequency domain
    # ----------------

    f_rel = fftfreq(dim_col, d=dt)
    f = f_rel + f_c
        
    # ideal differentiator
    H_id = 1j * DphiL_id(f_rel) / 2
    
    # MZI TF
    H = (1j * np.exp(-1j * DphiL(f) / 2) * np.sin(DphiL(f) / 2))


    matrix_der = np.zeros(img_array.shape) # matrix for the derivative at each pixel
    
    # ----------------
    # DOB (device operation bandwidth)
    # ----------------

    # it's the width of the band around the carrier in which the deviation of the ideal MZI from the ideal differentiator is within a certain tolerance

    f_id = f - df_c # it is defined for ideal f

    selection = (f_id > f_id[1]) # only positive frequencies and not the first ones (too close to 0)

    H_MZI_id = (1j * np.exp(-1j * DphiL_id(f_id) / 2) * np.sin(DphiL_id(f_id) / 2)) # it is defined for ideal f

    deviation = np.abs(np.abs(H_MZI_id[selection]) - np.abs(H_id[selection])) / np.abs(H_id[selection])

    tol_DOB = 0.01

    try:
        DOB = 2 * (f_id[selection][np.where(deviation <= tol_DOB)][-1] - f_c + df_c) # maximum frequency at which the deviation is smaller than the tolerance
    except IndexError: # in case it doesn't find it, the DOB is smaller than the frequency resolution
        DOB = 2 * (f_id[selection][0] - f_c + df_c) # upper bound
        print("Warning: DOB smaller than the frequency resolution, the value used is an upper bound")

    # ----------------
    # Some parameters
    # ----------------
    BWs = [] # it'll contain the bandwidth of the signal for each row
    ADs = [] # it'll contain the average deviation in the temporal output for each row

    # Quality factor
    Q = f_c/DOB

   
    # --------------------------------
    # Row elaboration
    # --------------------------------

    for i_row in range(dim_row):

        # the input waveform is the serialized stream of a row
        E_in = img_array[i_row] / np.sqrt(2) # half power since the carrier will create two copies (positive and negative side of the frequency axis)

        # ----------------
        # Fourier transform
        # ----------------

        E_in_fft = fft(E_in)
        E_out_fft = H * E_in_fft
        E_out = ifft(E_out_fft)
        
        E_in_der = ifft(1j * 2 * np.pi * f_rel * E_in_fft) # ideal derivative
        I_id = np.abs(E_in_der) ** 2

        I = np.abs(E_out) ** 2
        I_fft = fft(I)
        
        # ----------------
        # Maximum bandwidth and average deviations
        # ----------------
        
        BW_signal = bandwidth_threshold(fftshift(f), fftshift(np.abs(E_in_fft))) # both shifted, so that the peak is at the center of the axis
        BWs.append(BW_signal)
        
        # in a uniform row there are no edges, so the ideal derivative is zero and the average deviation is not defined
        is_uniform_row = np.ptp(E_in) == 0

        if not is_uniform_row:
            I_norm = I / np.max(I)
            I_id_norm = I_id / np.max(I_id)
            ADs.append(trapezoid(np.abs(I_norm - I_id_norm), x=t)/T) # calculation of the average deviation on the whole row

        # ----------------
        # Derivative
        # ----------------

        matrix_der[i_row] = I # for each pixel of the row the derivative is I_out = |E_out|^2

        # removing the edges of the images since at the starting and ending of the waveform the amplitude goes from 0 to x, so it sees an high derivative
        matrix_der[i_row, 0:3] = 0
        matrix_der[i_row, dim_col-3:dim_col] = 0

        # ----------------
        # Control plots
        # ----------------

        if is_control_plots and i_row == 1 and are_rows and not is_uniform_row: # change here which row to control

            E_in_fft = fftshift(E_in_fft)
            E_out_fft = fftshift(E_out_fft)

            fig, axs = plt.subplots(2, 2, figsize=(12, 8))

            # Temporal input
            axs[0, 0].plot(t, E_in, 
                           label='$E_\\mathrm{in}$', color='b')
            axs[0, 0].set_xlabel('$t$ [s]')
            axs[0, 0].set_ylabel('$E_\\mathrm{in}$')
            axs[0, 0].set_xlim(0, T)
            axs[0, 0].legend(fontsize=13)
            axs[0, 0].grid()
            axs[0, 0].set_title("Temporal input")

            # Temporal output
            axs[0, 1].plot(t, I_norm, 
                           label='$I_\\mathrm{out}$', color='b')
            axs[0, 1].plot(t, I_id_norm, 
                           label='$I_{out, id}$', color='g')
            axs[0, 1].plot(t, np.abs(I_norm - I_id_norm), 
                           label='$|I_\\mathrm{out} - I_{out, id}|$', color='k')
            axs[0, 1].set_xlabel('$t$ [s]')
            axs[0, 1].set_ylabel('$I_{norm}$')
            axs[0, 1].set_xlim(0, T)
            axs[0, 1].set_ylim(-1, 3)
            axs[0, 1].legend(fontsize=13)
            axs[0, 1].grid()
            axs[0, 1].set_title("Temporal output")

            # Input FFT
            plot_E_in_fft, = axs[1, 0].plot(fftshift(f), np.abs(E_in_fft), label="$|\\mathcal{F}\\{E_{in}(t)\\}|$", color='b')

            plot_H, = axs[1, 0].plot(fftshift(f), np.abs(fftshift(H)) * max(np.abs(E_in_fft)) / max(np.abs(H)), 
                label="$|H(f)|$", color='r')
            
            plot_H_id, = axs[1, 0].plot(fftshift(f), np.abs(fftshift(H_id)) * max(np.abs(E_in_fft)) / max(np.abs(H)), 
                label="$|H_{id}(f)|$", color='g')

            axs[1, 0].axvline(DOB/2+f_c, color='k', linestyle='-.', linewidth=1)
            line_DOB = axs[1, 0].axvline(-DOB/2+f_c, color='k', linestyle='-.', linewidth=1, label=f'DOB')

            axs[1, 0].set_xlabel('$f$ [Hz]')
            axs[1, 0].set_ylabel('Amplitude')
            axs[1, 0].set_xlim(f_c - f_s/10, f_c + f_s/10)
            axs[1, 0].set_ylim(-1, 1 + max(np.abs(E_in_fft)))

            axs[1, 0].axvline(f_c - BW_signal/2, color='cornflowerblue', linestyle='--', linewidth=2,)
            axs[1, 0].axvline(f_c + BW_signal/2,  color='cornflowerblue', linestyle='--', linewidth=2)
            area_BW_signal = axs[1, 0].axvspan(f_c - BW_signal/2, f_c + BW_signal/2, color='cornflowerblue', alpha=0.15, label='ISB')

            
            axs[1, 0].legend(handles=[plot_E_in_fft, plot_H, plot_H_id, line_DOB, area_BW_signal], fontsize=13)
            
            axs[1, 0].grid()
            axs[1, 0].set_title("Input amplitude spectrum")

            # Output FFT
            I_id_fft = fft(I_id)

            axs[1, 1].plot(fftshift(f), np.abs(fftshift(I_fft) / max(np.abs(fftshift(I_fft)))),
                label="$|\\mathcal{F}\\{I_\\mathrm{out}(t)\\}|$", color='b')

            axs[1, 1].plot(fftshift(f), np.abs(fftshift(I_id_fft) / max(np.abs(fftshift(I_id_fft)))),
                label="$|\\mathcal{F}\\{I_{out, id}(t)\\}|$", color='g')

            axs[1, 1].set_xlabel('$f$ [Hz]')
            axs[1, 1].set_ylabel('Amplitude')
            axs[1, 1].set_xlim(f_c - f_s, f_c + f_s)
            axs[1, 1].legend(fontsize=13)
            axs[1, 1].grid()
            axs[1, 1].set_title("Output amplitude spectrum")

            plt.tight_layout()
            plt.show()
            
    # ----------------
    # Fidelity (overlap integral)
    # ----------------
    
    BW_max = np.max(BWs)
    f_rel_dense = np.linspace(-BW_max/2, BW_max/2, int(1e4))
    f_dense = f_rel_dense + f_c

    # redefined TFs on a denser axis to accurately evaluate the integral
    H_id_dense = 1j * DphiL_id(f_rel_dense) / 2
    H_MZI_dense = (1j * np.exp(-1j * DphiL(f_dense) / 2) * np.sin(DphiL(f_dense) / 2))

    overlap = trapezoid(np.conj(H_id_dense) * H_MZI_dense, x=f_dense)

    norm_H_id = trapezoid(np.abs(H_id_dense)**2, x=f_dense)

    norm_H = trapezoid(np.abs(H_MZI_dense)**2, x=f_dense)

    F = np.abs(overlap)**2/(norm_H_id * norm_H)

    # ----------------
    # Final threshold to highlight the edges
    # ----------------
    
    mu = np.mean(matrix_der)
    sigma = np.std(matrix_der)

    thresh = mu + 1.7 * sigma # adaptive technique based on the average and std of the specific image
    
    matrix_der = matrix_der > thresh
    
    # ----------------
    # Figures of merit
    # ----------------
    
    BW_max_to_detun_min = 20 # empirical minimum of BW_max/df, below it the signal sees the detuning of the notch

    BW_max_to_detun = 1e10
    BW_max_to_detun_check = "OK"
    if df_c != 0:
        BW_max_to_detun = BW_max/df
        if BW_max_to_detun < BW_max_to_detun_min:
            BW_max_to_detun_check = "f_s too low"
    
    DOB_util_max = BW_max/DOB
    DOB_util_max_check = "OK"
    if DOB_util_max > 1:
        DOB_util_max_check = "f_s too high"


    # -------------------------------------------------------------------------
    # Prints
    # -------------------------------------------------------------------------

    if are_rows:
        print("===== Circuit and input parameters ========")
        print(f"DL = {DL_id:.3g} + {dDL:.3g} m")
        print(f"f_c = {f_c_id:.3g} + {df_c:.3g} Hz")
        print(f"f_s = {f_s:.3g} Hz")
        print(f"FSR = {FSR:.3g} Hz")
        print(f"Dm uncompensated = {Dm_n_uncomp:.7g} + {Dm_uncomp - Dm_n_uncomp:.7g} = {Dm_uncomp:.7g}")
        print(f"Dm compensated ={Dm_uncomp:.7g} + {Dm_comp-Dm_uncomp:.7g} = {Dm_comp:.7g}")
        print(f"dT used to compensate = {dT:.3g} °C")
        print(f"BW_mean, BW_std = {np.mean(BWs):.3g} Hz, {np.std(BWs):.3g} Hz")
        print(f"BW_min, BW_max = {np.min(BWs):.3g} Hz, {BW_max:.3g} Hz")

        print("\n===== Maximum errors ====================")
        print(f"Q = {Q:.3g}")
        print(f"dDL = {DL_id/Q:.3g} m")
        print(f"dn_eff = {n_eff_id/Q:g}")
        print(f"Waveguide width (w = 450nm) = {n_eff_id/(Q*dn_eff_dw):.3g} nm")
        print(f"Waveguide thickness (h = 220nm) = {n_eff_id/(Q*dn_eff_dh):.3g} nm")
        print(f"Temperature dT ={n_eff_id/(K_th*Q):.3g} °C")
        
        print("\n===== Figures of merit ====================")
        print(f"BW_max/df = {BW_max_to_detun:.3g} -> {BW_max_to_detun_check}")
        print(f"BW_max/DOB = {DOB_util_max:.3g} -> {DOB_util_max_check}")
        print(f"df/FSR = {df/FSR:.3g}")
        print(f"df/DOB = {df/DOB:.3g}")
        print(f"DOB = {DOB:.3g} Hz")
        
        print("\n===== Performance metrics ===================")
        print(f"F = {F:.3g}")
        print(f"AD_mean, AD_std = {100*np.mean(ADs):.3g}%, {100*np.std(ADs):.3g}%")
        print(f"AD_min, AD_max = {100*np.min(ADs):.3g}%, {100*np.max(ADs):.3g}%")
    
    return matrix_der




# =============================================================================
# EXECUTION
# =============================================================================

# the image is processed twice, the second time being the image transposed, since only rows are processed
# (only when the file is run as a script, so that other scripts can import its parameters and functions)

if __name__ == "__main__":

    final_matrix = np.where(
        compute_derivative(img_array, are_rows=True)[..., None], # actual rows
        red_color,
        original_image
    ) # final_matrix is the original image in which the red color is applied using as mask the variable matrix_der

    final_matrix = np.where(
        compute_derivative(img_array.transpose(), are_rows=False).transpose()[..., None], # rows of the transposed image, so columns of the original image
        red_color,
        final_matrix
    ) # same, but with the image already colored in red for the rows, so we add the processing for the columns


    if is_image_plot:
        plt.figure()
        plt.imshow(final_matrix)
        plt.axis('off')
        plt.show()



