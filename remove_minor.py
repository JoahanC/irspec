
from astropy.io import ascii
from astropy.table import Table 
from astropy.io import fits
from continuum_sub_test import trim_spec, cut_line, fit_continuum
from lmfit.models import GaussianModel
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import scipy.special as sp
from fitfuncs import * 
from plotparams import PlotParams
from emission_io import read_line_params

line_dict = read_line_params()
def gaussian(x, amplitude, center, sigma):
    s2pi = np.sqrt(2*np.pi)
    return (amplitude/(s2pi*sigma)) * np.exp(-(1.0*x-center)**2 / (2*sigma**2))

def gaussian_integral(amplitude, center, sigma, x_1, x_2):
    return amplitude * (sp.erf((center - x_1)) / (np.sqrt(2) * sigma) - sp.erf((center - x_2)) / (np.sqrt(2) * sigma)) / 2

def bkg_integral(m, b, x_1, x_2):
    return m * (x_2 ** 2 / 2 - x_1 ** 2 / 2) + b * (x_2 - x_1)

def evaluate_spaxel_fit(name, redshift, test_pixel):
    
    i = test_pixel[1]
    j = test_pixel[0]
    # Define dictionary of emission line values
    line_dict = {"[NeV]": [4, "long", [24.1, 24.5], 24.316],
                "[NeV]_14": [3, "medium", [14.28, 14.36], [14.23, 14.40], 14.3217, [0.03, 0.1, 0.0125, 0.02]],
                "[H_2_S_1]": [3, "long", [16.99, 17.06], [16.98, 17.07], 17.03, [0.03, 0.1, 0.0125, 0.02]], 
                "[H_2_S_3]": [2, "medium", [9.63, 9.71], [9.57, 9.8], 9.6654, [0.03, 0.1, 0.0125, 0.02]], 
                "[OIV]": [4, "long", [25.62, 26.08], 25.88], 
                "MgV": [1, "medium", [5.57, 5.656], 5.609],
                "[NeIII]": [3, "long", [15.5251, 15.6051], [15.38, 15.67], 15.5551, [0.03, 0.001, 0.005, 0.01]],
                "[SIV]": [2, "long", [10.3, 10.63], [10.00, 10.87], 10.509, [0.03, 0.1, 0.0125, 0.02]],
                "[NeII]": [3, "short", [12.77, 12.84], [12.73, 12.89], 12.813550, [0.03, 0.1, 0.0125, 0.02]],
                "[FeII]": [1, "short", [5.29, 5.37], [5.25, 5.43], 5.3396, [0.01, 0.5, 0.0125, 0.02]]}
    
    # Access fit data
    fitdata = ascii.read(f"./../diagnostic_plots/dynamic_multicomponent/{name}/AMP_fit.dat", format="ipac")  

    # Access datacube and assemble wavelength and flux values
    datacube = f"./../input_data/IR23128-S/ch{line_dict[name][0]}-{line_dict[name][1]}_s3d.fits"
    hdul = fits.open(datacube)
    data = hdul[1].data
    wave0 = hdul[1].header["CRVAL3"]
    dwave = hdul[1].header["CDELT3"]
    steps = np.arange(len(data))
    deltawave = steps * dwave
    wavelengths = (wave0 + deltawave) / (1 + redshift)

    # Extract and process spectrum for spaxel
    full_spec_1 = data[:, i, j]
    trim_wavelengths, trim_flux = trim_spec(wavelengths, full_spec_1, line_dict[name][3][0], line_dict[name][3][1])
    new_wave, new_flux, transition_idx = cut_line(trim_wavelengths, trim_flux, line_dict[name][2][0], line_dict[name][2][1])
    popt, pcov = fit_continuum(new_wave, new_flux)
    bkg_sub_fluxes = trim_flux - OneDPolynomial(trim_wavelengths, popt[0], popt[1])

    # Get fit data
    for idx, _ in enumerate(fitdata["XPIX"]):
        if fitdata["XPIX"][idx] == i:
            if fitdata["YPIX"][idx] == j:
                g1_amp = fitdata["G1AMP"][idx]
                g1_cen = fitdata["G1CEN"][idx]
                g1_sigma = fitdata["G1SIGMA"][idx]
                g2_amp = fitdata["G2AMP"][idx]
                g2_cen = fitdata["G2CEN"][idx]
                g2_sigma = fitdata["G2SIGMA"][idx]
                g3_amp = fitdata["G3AMP"][idx]
                g3_cen = fitdata["G3CEN"][idx]
                g3_sigma = fitdata["G3SIGMA"][idx]

    # Reset plot parameters
    pltparams = PlotParams(palatte="dark", scaling="paper")



    # Initialize plot
    fig, ax = plt.subplots()

    # Find peaks and label them
    peaks, _ = find_peaks(bkg_sub_fluxes, height=np.max(bkg_sub_fluxes)*0.2)
    trim_peaks = []
    for peak in peaks:
        if trim_wavelengths[peak] > line_dict[name][2][0] and trim_wavelengths[peak] < line_dict[name][2][1]:
            trim_peaks.append(peak)

    # Plot spectrum
    ax.plot(trim_wavelengths, bkg_sub_fluxes, c="white", alpha=0.75)

    # Iterate through gaussian components and evaluate each
    if not np.isnan(g1_amp):
        
        # Initialize and plot gaussian 
        g1 = GaussianModel()
        params_g1 = g1.guess(bkg_sub_fluxes, x=trim_wavelengths)
        params_g1.update(g1.make_params(center=dict(value=g1_cen),
                                                sigma=dict(value=g1_sigma),
                                                amplitude=dict(value=g1_amp)))
        g1_predict = g1.eval(params_g1, x=trim_wavelengths)
        g_predict = g1_predict
        ax.plot(trim_wavelengths, g1_predict, c="cyan", label=f"Gaussian 1", lw=2)
        
        # Estimate SNR calculation region
        g1_fwhm = g1_sigma * 2*np.sqrt(2*np.log(2))
        snr_wavelengths = np.linspace(g1_cen - g1_fwhm, g1_cen + g1_fwhm, 1000)
        g1_flux = g1.eval(params_g1, x=snr_wavelengths)
        bkg_flux = OneDPolynomial(snr_wavelengths, popt[0], popt[1])
        
        # Estimate SNR
        line_signal = np.trapz(g1_flux + bkg_flux, snr_wavelengths)
        bkg_signal = np.trapz(bkg_flux, snr_wavelengths)
        snr = (line_signal + bkg_signal) / np.sqrt(bkg_signal)
        
        # Read out SNR
        print(f"Gaussian Component 1: {round(g1_cen, 3)}")
        print(f"SNR: {snr}")
        if snr < 5:
            print("Rejecting component.")
        else:
            print("Accepting component.")

    if not np.isnan(g2_amp):
        g2 = GaussianModel()
        params_g2 = g2.guess(bkg_sub_fluxes, x=trim_wavelengths)
        params_g2.update(g2.make_params(center=dict(value=g2_cen),
                                                sigma=dict(value=g2_sigma),
                                                amplitude=dict(value=g2_amp)))
        g2_predict = g2.eval(params_g2, x=trim_wavelengths)
        g_predict += g2_predict
        ax.plot(trim_wavelengths, g2_predict, c="lime", label=f"Gaussian 2", lw=0.9)
        
        # Estimate SNR calculation region
        g2_fwhm = g2_sigma * 2*np.sqrt(2*np.log(2))
        snr_wavelengths = np.linspace(g2_cen - g2_fwhm, g2_cen + g2_fwhm, 1000)
        g2_flux = g2.eval(params_g2, x=snr_wavelengths)
        bkg_flux = OneDPolynomial(snr_wavelengths, popt[0], popt[1])
        
        # Estimate SNR
        line_signal = np.trapz(g2_flux + bkg_flux, snr_wavelengths)
        bkg_signal = np.trapz(bkg_flux, snr_wavelengths)
        snr = (line_signal + bkg_signal) / np.sqrt(bkg_signal)
        
        print(f"Gaussian Component 2: {round(g2_cen, 3)}")
        print(f"SNR: {snr}")
        if snr < 5:
            print("Rejecting component.")
        else:
            print("Accepting component.")

    if not np.isnan(g3_amp):
        g3 = GaussianModel()
        params_g3 = g3.guess(bkg_sub_fluxes, x=trim_wavelengths)
        params_g3.update(g3.make_params(center=dict(value=g3_cen),
                                                sigma=dict(value=g3_sigma),
                                                amplitude=dict(value=g3_amp)))
        g3_predict = g3.eval(params_g3, x=trim_wavelengths)
        g_predict += g3_predict
        ax.plot(trim_wavelengths, g3_predict, c="orangered", label=f"Gaussian 3", lw=0.9)
        
        # Estimate SNR calculation region
        g3_fwhm = g3_sigma * 2*np.sqrt(2*np.log(2))
        snr_wavelengths = np.linspace(g3_cen - g3_fwhm, g3_cen + g3_fwhm, 1000)
        g3_flux = g3.eval(params_g3, x=snr_wavelengths)
        bkg_flux = OneDPolynomial(snr_wavelengths, popt[0], popt[1])
        
        # Estimate SNR
        line_signal = np.trapz(g3_flux + bkg_flux, snr_wavelengths)
        bkg_signal = np.trapz(bkg_flux, snr_wavelengths)
        snr = (line_signal + bkg_signal) / np.sqrt(bkg_signal)
        
        print(f"Gaussian Component 3: {round(g3_cen, 3)}")
        print(f"SNR: {snr}")
        if snr < 5:
            print("Rejecting component.")
        else:
            print("Accepting component.")

    # Plot cumulative fit
    ax.plot(trim_wavelengths, g_predict, c="gold", ls="dashed", label=f"All Gaussians", lw=0.9)
    ax.legend()
    plt.savefig(f"./../diagnostic_plots/dynamic_multicomponent/{name}/test_spaxels/{i}_{j}_trimtest.pdf", bbox_inches="tight")
    plt.close()

def evaluate_gaussian_significance(m, b, amplitude, center, sigma):
    g = GaussianModel()
    params_g = g.make_params(center=dict(value=0),
                                sigma=dict(value=0),
                                amplitude=dict(value=0))
    params_g.update(g.make_params(center=dict(value=center),
                                    sigma=dict(value=sigma),
                                    amplitude=dict(value=amplitude)))
    
    fwhm = sigma * 2*np.sqrt(2*np.log(2))
    snr_wavelengths = np.linspace(center - fwhm, center + fwhm, 1000)
    g_flux = g.eval(params_g, x=snr_wavelengths)
    bkg_flux = OneDPolynomial(snr_wavelengths, m, b)
    
    # Estimate SNR
    line_signal = np.trapz(g_flux + bkg_flux, snr_wavelengths)
    bkg_signal = np.trapz(bkg_flux, snr_wavelengths)
    snr = (line_signal + bkg_signal) / np.sqrt(bkg_signal)
    return snr

def remove_insignificant_gaussians(name, snr_cutoff):
    
    line_dict = {"[NeV]": [4, "long", [24.1, 24.5], 24.316],
                "[NeV]_14": [3, "medium", [14.28, 14.36], [14.23, 14.40], 14.3217, [0.03, 0.1, 0.0125, 0.02]],
                "[H_2_S_1]": [3, "long", [16.99, 17.06], [16.98, 17.07], 17.03, [0.03, 0.1, 0.0125, 0.02]], 
                "[H_2_S_3]": [2, "medium", [9.63, 9.71], [9.57, 9.8], 9.6654, [0.03, 0.1, 0.0125, 0.02]], 
                "[OIV]": [4, "long", [25.62, 26.08], 25.88], 
                "MgV": [1, "medium", [5.57, 5.656], 5.609],
                "[NeIII]": [3, "long", [15.5251, 15.6051], [15.38, 15.67], 15.5551, [0.03, 0.001, 0.005, 0.01]],
                "[SIV]": [2, "long", [10.3, 10.63], [10.00, 10.87], 10.509, [0.03, 0.1, 0.0125, 0.02]],
                "[NeII]": [3, "short", [12.77, 12.84], [12.73, 12.89], 12.813550, [0.03, 0.1, 0.0125, 0.02]],
                "[FeII]": [1, "short", [5.29, 5.37], [5.25, 5.43], 5.3396, [0.01, 0.5, 0.0125, 0.02]]}
    
    g1_amps = []
    g2_amps = []
    g3_amps = []
    g1_vels = []
    g2_vels = []
    g3_vels = []
    g1_sigmas = []
    g2_sigmas = []
    g3_sigmas = []
    
    # Access fit data
    datacube = f"./../input_data/IR23128-S/ch{line_dict[name][0]}-{line_dict[name][1]}_s3d.fits"
    hdul = fits.open(datacube)
    data = hdul[1].data
    wave0 = hdul[1].header["CRVAL3"]
    dwave = hdul[1].header["CDELT3"]
    steps = np.arange(len(data))
    deltawave = steps * dwave
    wavelengths = (wave0 + deltawave) / (1 + redshift)
    fitdata = ascii.read(f"./../diagnostic_plots/dynamic_multicomponent/{name}/fit.dat", format="ipac")
    
    # Loop through all gaussians
    for idx, _ in enumerate(fitdata["XPIX"]):
        
        
        # Check is gauss1 is nan
        if not np.isnan(fitdata["G1AMP"][idx]):
            
            # Extract and process spectrum for spaxel
            full_spec_1 = data[:, fitdata["XPIX"][idx], fitdata["YPIX"][idx]]
            trim_wavelengths, trim_flux = trim_spec(wavelengths, full_spec_1, line_dict[name][3][0], line_dict[name][3][1])
            new_wave, new_flux, transition_idx = cut_line(trim_wavelengths, trim_flux, line_dict[name][2][0], line_dict[name][2][1])
            popt, pcov = fit_continuum(new_wave, new_flux)
            
            # Check if gauss2 is nan
            if not np.isnan(fitdata["G2AMP"][idx]):
                # Check if gauss3 is nan
                if not np.isnan(fitdata["G3AMP"][idx]):
                    # Check all gauss significance
                    snr1 = evaluate_gaussian_significance(popt[0], popt[1], fitdata["G1AMP"][idx], fitdata["G1CEN"][idx], fitdata["G1SIGMA"][idx])
                    snr2 = evaluate_gaussian_significance(popt[0], popt[1], fitdata["G2AMP"][idx], fitdata["G2CEN"][idx], fitdata["G2SIGMA"][idx])
                    snr3 = evaluate_gaussian_significance(popt[0], popt[1], fitdata["G3AMP"][idx], fitdata["G3CEN"][idx], fitdata["G3SIGMA"][idx])
                    
                    valid_gaussian_idicies = []
                    if snr1 >= snr_cutoff:
                        valid_gaussian_idicies.append("1")
                    if snr2 >= snr_cutoff:
                        valid_gaussian_idicies.append("2")
                    if snr3 >= snr_cutoff:
                        valid_gaussian_idicies.append("3")
                    if len(valid_gaussian_idicies) == 3:
                        g1_amps.append(fitdata["G1AMP"][idx])
                        g2_amps.append(fitdata["G2AMP"][idx])
                        g3_amps.append(fitdata["G3AMP"][idx])
                        g1_vels.append(fitdata["G1CEN"][idx])
                        g2_vels.append(fitdata["G2CEN"][idx])
                        g3_vels.append(fitdata["G3CEN"][idx])
                        g1_sigmas.append(fitdata["G1SIGMA"][idx])
                        g2_sigmas.append(fitdata["G2SIGMA"][idx])
                        g3_sigmas.append(fitdata["G3SIGMA"][idx])
                    if len(valid_gaussian_idicies) == 2:
                        g1_amps.append(fitdata[f"G{valid_gaussian_idicies[0]}AMP"][idx])
                        g2_amps.append(fitdata[f"G{valid_gaussian_idicies[1]}AMP"][idx])
                        g3_amps.append(np.nan)
                        g1_vels.append(fitdata[f"G{valid_gaussian_idicies[0]}CEN"][idx])
                        g2_vels.append(fitdata[f"G{valid_gaussian_idicies[1]}CEN"][idx])
                        g3_vels.append(np.nan)
                        g1_sigmas.append(fitdata[f"G{valid_gaussian_idicies[0]}SIGMA"][idx])
                        g2_sigmas.append(fitdata[f"G{valid_gaussian_idicies[1]}SIGMA"][idx])
                        g3_sigmas.append(np.nan)
                    if len(valid_gaussian_idicies) == 1:
                        g1_amps.append(fitdata[f"G{valid_gaussian_idicies[0]}AMP"][idx])
                        g2_amps.append(np.nan)
                        g3_amps.append(np.nan)
                        g1_vels.append(fitdata[f"G{valid_gaussian_idicies[0]}CEN"][idx])
                        g2_vels.append(np.nan)
                        g3_vels.append(np.nan)
                        g1_sigmas.append(fitdata[f"G{valid_gaussian_idicies[0]}SIGMA"][idx])
                        g2_sigmas.append(np.nan)
                        g3_sigmas.append(np.nan)
                    if len(valid_gaussian_idicies) == 0:
                        g1_amps.append(np.nan)
                        g2_amps.append(np.nan)
                        g3_amps.append(np.nan)
                        g1_vels.append(np.nan)
                        g2_vels.append(np.nan)
                        g3_vels.append(np.nan)
                        g1_sigmas.append(np.nan)
                        g2_sigmas.append(np.nan)
                        g3_sigmas.append(np.nan)
                    #print(snr1, snr2, snr3)
                    
                # Check gauss1 and gauss2 significance
                else:
                    snr1 = evaluate_gaussian_significance(popt[0], popt[1], fitdata["G1AMP"][idx], fitdata["G1CEN"][idx], fitdata["G1SIGMA"][idx])
                    snr2 = evaluate_gaussian_significance(popt[0], popt[1], fitdata["G2AMP"][idx], fitdata["G2CEN"][idx], fitdata["G2SIGMA"][idx])
                    #print(snr1, snr2)
                    valid_gaussian_idicies = []
                    if snr1 >= snr_cutoff:
                        valid_gaussian_idicies.append("1")
                    if snr2 >= snr_cutoff:
                        valid_gaussian_idicies.append("2")
                    
                    if len(valid_gaussian_idicies) == 2:
                        g1_amps.append(fitdata[f"G{valid_gaussian_idicies[0]}AMP"][idx])
                        g2_amps.append(fitdata[f"G{valid_gaussian_idicies[1]}AMP"][idx])
                        g1_vels.append(fitdata[f"G{valid_gaussian_idicies[0]}CEN"][idx])
                        g2_vels.append(fitdata[f"G{valid_gaussian_idicies[1]}CEN"][idx])
                        g1_sigmas.append(fitdata[f"G{valid_gaussian_idicies[0]}SIGMA"][idx])
                        g2_sigmas.append(fitdata[f"G{valid_gaussian_idicies[1]}SIGMA"][idx])
                    if len(valid_gaussian_idicies) == 1:
                        g1_amps.append(fitdata[f"G{valid_gaussian_idicies[0]}AMP"][idx])
                        g2_amps.append(np.nan)
                        g1_vels.append(fitdata[f"G{valid_gaussian_idicies[0]}CEN"][idx])
                        g2_vels.append(np.nan)
                        g1_sigmas.append(fitdata[f"G{valid_gaussian_idicies[0]}SIGMA"][idx])
                        g2_sigmas.append(np.nan)
                    if len(valid_gaussian_idicies) == 0:
                        g1_amps.append(np.nan)
                        g2_amps.append(np.nan)
                        g1_vels.append(np.nan)
                        g2_vels.append(np.nan)
                        g1_sigmas.append(np.nan)
                        g2_sigmas.append(np.nan)
                    
                    g3_amps.append(fitdata["G3AMP"][idx])
                    g3_vels.append(fitdata["G3CEN"][idx])
                    g3_sigmas.append(fitdata["G3SIGMA"][idx])
                    
            # Only check gauss1 significance
            else:
                snr1 = evaluate_gaussian_significance(popt[0], popt[1], fitdata["G1AMP"][idx], fitdata["G1CEN"][idx], fitdata["G1SIGMA"][idx])
                if np.isnan(snr1):
                    g1_amps.append(np.nan)
                    g1_vels.append(np.nan)
                    g1_sigmas.append(np.nan)
                if snr1 < snr_cutoff:
                    g1_amps.append(np.nan)
                    g1_vels.append(np.nan)
                    g1_sigmas.append(np.nan)
                if snr1 >= snr_cutoff:
                    g1_amps.append(fitdata["G1AMP"][idx])
                    g1_vels.append(fitdata["G1CEN"][idx])
                    g1_sigmas.append(fitdata["G1SIGMA"][idx])
                
                g2_amps.append(fitdata["G2AMP"][idx])
                g3_amps.append(fitdata["G3AMP"][idx])
                g2_vels.append(fitdata["G2CEN"][idx])
                g3_vels.append(fitdata["G3CEN"][idx])
                g2_sigmas.append(fitdata["G2SIGMA"][idx])
                g3_sigmas.append(fitdata["G3SIGMA"][idx])
        # All nans, no action necessary
        else:
            g1_amps.append(fitdata["G1AMP"][idx])
            g2_amps.append(fitdata["G2AMP"][idx])
            g3_amps.append(fitdata["G3AMP"][idx])
            g1_vels.append(fitdata["G1CEN"][idx])
            g2_vels.append(fitdata["G2CEN"][idx])
            g3_vels.append(fitdata["G3CEN"][idx])
            g1_sigmas.append(fitdata["G1SIGMA"][idx])
            g2_sigmas.append(fitdata["G2SIGMA"][idx])
            g3_sigmas.append(fitdata["G3SIGMA"][idx])
    
    #print(len(g1_amps), len(g2_amps), len(g3_amps), len(g1_vels), len(g2_vels), len(g3_vels), len(g1_sigmas), len(g2_sigmas), len(g3_sigmas))
    best_fit_values = [fitdata["XPIX"], fitdata["YPIX"], g1_amps, g1_vels, g1_sigmas, g2_amps, g2_vels, g2_sigmas, g3_amps, g3_vels, g3_sigmas]
    fitparams = Table(best_fit_values, names=("XPIX", "YPIX", "G1AMP", "G1CEN", "G1SIGMA", "G2AMP", "G2CEN", "G2SIGMA", "G3AMP", "G3CEN", "G3SIGMA"))
    fitparams.write(f"./../diagnostic_plots/dynamic_multicomponent/{name}/SNRCUT_fit.dat", format="ipac", overwrite=True)

# Define emission line parameters and spaxel of interest
name = "[NeIII]"
redshift = 0.044601
test_pixel = (22, 30) #(Y, X))
snr_cutoff = 2

#evaluate_spaxel_fit(name, redshift, test_pixel)
remove_insignificant_gaussians(name, snr_cutoff)