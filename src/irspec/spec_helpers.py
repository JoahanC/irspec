import numpy as np
from scipy.optimize import curve_fit

from irspec.fitfuncs import *

def trim_spec(wavelengths, fluxes, low_wave, high_wave):
    """ 
    Truncates a (wavelength, flux) pair of arrays to within
    two wavelength values.
    
    Arguments
    ---------
    wavelengths : array
        A 2D array of n elements corresponding to wavelength values
    fluxes : array 
        A 2D array of n elements corresponding to wavelength values
    low_wave : float
        The lower wavelength value
    high_wave : float 
        The upper wavelength value
    """
    min_diff_array = np.absolute(wavelengths-low_wave)
    max_diff_array = np.absolute(wavelengths-high_wave)
    low_idx = min_diff_array.argmin()
    high_idx = max_diff_array.argmin() 
    new_wave = wavelengths[low_idx:high_idx]
    new_flux = fluxes[low_idx:high_idx]
    return new_wave, new_flux


def cut_line(wavelengths, fluxes, low_wave, high_wave):
    """
    Cuts out a segment of a (wavelength, flux) pair of arrays 
    and concatenates the remaining sections.
    
    Arguments
    ---------
    wavelengths : array 
        A 2D array of n elements corresponding to wavelength values
    fluxes : array 
        A 2D array of n elements corresponding to wavelength values
    low_wave : float
        The lower wavelength value
    high_wave : float
        The upper wavelength value
    """
    min_dif_array = np.absolute(wavelengths-low_wave)
    max_dif_array = np.absolute(wavelengths-high_wave)
    low_idx = min_dif_array.argmin()
    high_idx = max_dif_array.argmin() 
    new_wave = np.concatenate((wavelengths[0:low_idx], wavelengths[high_idx:-1]))
    new_flux = np.concatenate((fluxes[0:low_idx], fluxes[high_idx:-1]))
    return new_wave, new_flux, low_idx


def fit_continuum(wavelengths, fluxes):
    """ 
    Fits a linear continuum to a (wavelength, flux) pair of arrays.
    
    Arguments
    ---------
    wavelengths : array 
        A 2D array of n elements corresponding to wavelength values
    fluxes : array 
        A 2D array of n elements corresponding to wavelength values
    """
    
    # Define fitting function
    def continuum_function(wave, *pars):
        continuum = OneDPolynomial(wave, pars[0], pars[1])
        fluxes = continuum 
        return fluxes
    
    guess = [0,np.min(fluxes)]
    popt, pcov = curve_fit(continuum_function, wavelengths, fluxes, guess, maxfev=5000)
    return popt, pcov


def find_nans(fluxes):
    """ 
    Finds continuous gaps in a spectra that are population by 
    `np.nan` values.
    """
    indicies = []
    in_nan_field = False
    for idx, val in enumerate(fluxes):
        if np.isnan(val):
            if in_nan_field == False:
                start_idx = idx 
                in_nan_field = True
            if idx == len(fluxes) - 1:
                if in_nan_field == True:
                    indicies.append((start_idx, idx))
        else:
            if in_nan_field == True:
                indicies.append((start_idx, idx - 1))
                in_nan_field = False
                start_idx = None
    return indicies