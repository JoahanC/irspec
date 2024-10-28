import numpy as np
import matplotlib.pyplot as plt 
import astropy.io.fits as fits
from scipy.optimize import curve_fit

from fitfuncs import *



"""core = "S"
name = "[NeII]"
type = "triple"
multicomponent = True

line_dict = {"[NeV]": [4, "long", [24.1, 24.5], 24.316],
             "[NeV]_14": [3, "medium", [14.28, 14.35], 14.3217],
             "[H_2_S_3]": [2, "medium", [8.975, 9.01], [8.8, 9.2], 8.98], 
             "[OIV]": [4, "long", [25.62, 26.08], 25.88], 
             "MgV": [1, "medium", [5.57, 5.656], 5.609],
             "[NeIII]": [3, "long", [15.38, 15.67], 15.5551],
             "[NeII]": [3, "short", [12.77, 12.84], [12.73, 12.89], 12.81],
             "[FeII]": [1, "short", [5.30, 5.41], 5.3396]}

redshift = 0.044601
if core == "N":
    datacube = f"./../input_data/IR23128-N/ch{line_dict[name][0]}-{line_dict[name][1]}_s3d.fits"
if core == "S":
    datacube = f"./../input_data/IR23128-S/ch{line_dict[name][0]}-{line_dict[name][1]}_s3d.fits"""

def trim_spec(wavelengths, fluxes, line_low, line_high):
    min_dif_array = np.absolute(wavelengths-line_low)
    max_dif_array = np.absolute(wavelengths-line_high)
    low_idx = min_dif_array.argmin()
    high_idx = max_dif_array.argmin() 
    new_wave = wavelengths[low_idx:high_idx]
    new_flux = fluxes[low_idx:high_idx]
    return new_wave, new_flux


def cut_line(wavelengths, fluxes, line_low, line_high):
    min_dif_array = np.absolute(wavelengths-line_low)
    max_dif_array = np.absolute(wavelengths-line_high)
    low_idx = min_dif_array.argmin()
    high_idx = max_dif_array.argmin() 
    new_wave = np.concatenate((wavelengths[0:low_idx], wavelengths[high_idx:-1]))
    new_flux = np.concatenate((fluxes[0:low_idx], fluxes[high_idx:-1]))
    return new_wave, new_flux, low_idx


def fit_continuum(wavelengths, fluxes):
    
    def continuum_function(wave, *pars):
        continuum = OneDPolynomial(wave, pars[0], pars[1])
        fluxes = continuum 
        return fluxes
    guess = [0,np.min(fluxes)]
    popt, pcov = curve_fit(continuum_function, wavelengths, fluxes, guess, maxfev=5000)

    return popt, pcov


"""hdul = fits.open(datacube)
data = hdul[1].data
wave0 = hdul[1].header["CRVAL3"]
dwave = hdul[1].header["CDELT3"]

steps = np.arange(len(data))
deltawave = steps * dwave
wavelengths = (wave0 + deltawave) / (1 + redshift)

data_quality = hdul[3].data[0]
fluxes = np.zeros((np.shape(data)[1], np.shape(data)[2]))

fig, ax = plt.subplots()
for i in range(np.shape(data)[1]):
    for j in range(np.shape(data)[2]):
        
        full_spec_1 = data[:, i, j]
        if np.isnan(full_spec_1).any():
            print("NAN IN SPECTRUM")
            continue
        trim_wavelengths, trim_flux = trim_spec(wavelengths, full_spec_1, line_dict[name][3][0], line_dict[name][3][1])
        new_wave, new_flux = cut_line(trim_wavelengths, trim_flux, line_dict[name][2][0], line_dict[name][2][1])
        popt, pcov = fit_continuum(new_wave, new_flux)
        bkg_sub_fluxes = trim_flux - OneDPolynomial(trim_wavelengths, popt[0], popt[1])
        ax.plot(trim_wavelengths, bkg_sub_fluxes)
plt.show()"""