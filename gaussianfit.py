import numpy as np
import matplotlib.pyplot as plt
from cubespec import CubeSpec
from scipy.optimize import curve_fit


spec_obj = CubeSpec("./../", "param_files", "IR23128-S_0_single_param.txt", "input_data/IR23128-S/", redshift=0.044601, fit_dirname="IR23128S_AGN0", mode="AGN")
data = spec_obj.recall_data()

line_center = 15.556
x_low = line_center - 0.2
x_high = line_center + 0.2

min_dif_array = np.absolute(data["wave"]-x_low)
max_dif_array = np.absolute(data["wave"]-x_high)
low_idx = min_dif_array.argmin()
high_idx = max_dif_array.argmin() 

data_cutout = {"wave": data["wave"][low_idx:high_idx], "flux": data["flux"][low_idx:high_idx]}


def gaussian(x, A, x0, sig):
    return A*np.exp(-(x-x0)**2/(2*sig**2))

def continuum(x, b, c):
    return b*x+c

def multi_gaussian(x, *pars):
    g1 = gaussian(x, pars[0], pars[1], pars[2])
    #g2 = gaussian(x, pars[3], pars[4], pars[5])
    #g3 = gaussian(x, pars[6], pars[7], pars[8])
    cont = continuum(x, pars[3], pars[4])
    return g1 + cont


def calc_sigma(wavelength, flux):
    peak_idx = np.argmax(flux)
    pre_flux = flux[:peak_idx]
    post_flux = flux[peak_idx:]
    root2val = np.max(flux) / np.sqrt(2)
    pre_dif_arr = np.absolute(pre_flux - root2val)
    post_diff_arr = np.absolute(post_flux - root2val)
    pre_dif_idx = np.argmin(pre_dif_arr)
    post_dif_idx = np.argmin(post_diff_arr) + peak_idx
    return wavelength[post_dif_idx] - wavelength[pre_dif_idx]

def fit_cutout(wavelength, flux, line_center):

    guess = [np.nan, np.nan, np.nan, np.nan, np.nan]
    narrow_cut_low_idx = np.argmin(np.absolute(wavelength - (line_center - 0.01)))
    narrow_cut_high_idx = np.argmin(np.absolute(wavelength - (line_center + 0.01)))
    guess[0] = np.max(flux[narrow_cut_low_idx:narrow_cut_high_idx])
    guess[4] = np.average(flux)
    guess[1] = line_center
    guess[2] = 0.02 #calc_sigma(wavelength, flux)
    guess[3] = (flux[-1] - flux[0]) / np.average(wavelength)
    print(guess)
    popt, pcov = curve_fit(multi_gaussian, data_cutout["wave"], data_cutout["flux"], guess)#, 
                           #bounds = ((guess[0] * 0.1, line_center - 0.01, 0.005, -np.inf, -np.inf), (guess[0], line_center + 0.01, 0.03, np.inf, np.inf)))
    return popt, pcov

popt, pcov = fit_cutout(data_cutout["wave"], data_cutout["flux"], line_center)
residuals = (multi_gaussian(data_cutout["wave"], *popt) - data_cutout["flux"]) / data_cutout["flux"] * 100

plt.style.use('dark_background')
fig, ((ax1), (ax2)) = plt.subplots(nrows=2,ncols=1)
fig.subplots_adjust(hspace=0.0)
ax1.scatter(data_cutout["wave"], data_cutout["flux"], s=5)
ax1.plot(data_cutout["wave"], multi_gaussian(data_cutout["wave"], *popt), 'r--', linewidth=2, label='Fit')
ax2.plot(data_cutout["wave"], residuals)
#ax.set_xscale("log")
ax2.set_xlabel(r"$\lambda$ (microns)")
ax1.set_ylabel(r"$f_{\nu} (Jy)$")
ax2.set_ylabel("Residuals (%)")
std = np.nanstd(residuals) 
ax2.set_ylim(-4*std, 4*std)
ax1.set_yscale("log")
#ax2.set_yscale("log")
ax1.legend()
plt.show()