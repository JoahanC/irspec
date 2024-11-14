""" 
This scripts implements a dynamic multicomponent fitting approach.
"""
import numpy as np 
from continuum_sub_test import trim_spec, cut_line, fit_continuum
from scipy.signal import find_peaks


from fitfuncs import * 
from lmfit.models import GaussianModel, LinearModel
from astropy.table import Table
from tqdm import tqdm

import matplotlib.pyplot as plt
from astropy.io import fits

from plotparams import PlotParams
from emission_io import read_line_params
pltparams = PlotParams(scaling="presentation")



def main_routine(core, name, redshift, test_pixel=False, experimental_bkg=False):
    """ 
    The main fitting routine for multicomponent fitting
    """
    
    # Initialize dictionary of line values
    # LINENAME: [CH, SUBCH, [TRIMWINDOW], [FITWINDOW], LINECENTER, [DELTALINE, DELTAC1, FWHM1, FWHM2, FWHM3]]
    """line_dict = {"[NeV]": [4, "long", [24.1, 24.5], 24.316],
             "[NeV]_14": [3, "medium", [14.28, 14.36], [14.23, 14.40], 14.3217, [0.03, 0.1, 0.0125, 0.02]],
             "[H_2_S_1]": [3, "long", [16.99, 17.06], [16.98, 17.07], 17.03, [0.03, 0.1, 0.0125, 0.02]], 
             "[H_2_S_3]": [2, "medium", [9.63, 9.71], [9.57, 9.8], 9.6654, [0.03, 0.1, 0.0125, 0.02]], 
             "[OIV]": [4, "long", [25.62, 26.08], 25.88], 
             "MgV": [1, "medium", [5.57, 5.656], 5.609],
             "[NeIII]": [3, "long", [15.5251, 15.6051], [15.38, 15.67], 15.5551, [0.03, 0.001, 0.005, 0.01]],
             "[SIV]": [2, "long", [10.3, 10.63], [10.00, 10.87], 10.509, [0.03, 0.1, 0.0125, 0.02]],
             "[NeII]": [3, "short", [12.77, 12.84], [12.73, 12.89], 12.813550, [0.03, 0.1, 0.0125, 0.02]],
             "[FeII]": [1, "short", [5.29, 5.37], [5.25, 5.43], 5.3396, [0.03, 0.1, 0.001, 0.005, 0.01]]}"""
    line_dict = read_line_params()
    
    # Initialize file path variable
    if core == "N":
        datacube = f"./../input_data/IR23128-N/ch{line_dict[name][0]}-{line_dict[name][1]}_s3d.fits"
    if core == "S":
        datacube = f"./../input_data/IR23128-S/ch{line_dict[name][0]}-{line_dict[name][1]}_s3d.fits"
    
    # Access datacube and assemble wavelength and flux values
    hdul = fits.open(datacube)
    data = hdul[1].data
    wave0 = hdul[1].header["CRVAL3"]
    dwave = hdul[1].header["CDELT3"]

    steps = np.arange(len(data))
    deltawave = steps * dwave
    wavelengths = (wave0 + deltawave) / (1 + redshift)
    
    data_quality = hdul[3].data[0]
    fluxes = np.zeros((np.shape(data)[1], np.shape(data)[2]))
    base_array = np.zeros((np.shape(data)[1], np.shape(data)[2]))
    base_array2 = np.zeros((np.shape(data)[1], np.shape(data)[2]))
    
    # Test pixel indicies
    test_pixels = [(32, 10), (26,16), (31,19), (26,36), (25,21), (28,38), (25,22), (35,12), (29, 10), (29, 32), (22, 26), (26,22), (19, 26), (26, 19), (31, 14)]
    #test_pixels = [(26, 19)]
    #test_pixels = [(31, 14)]
    
    # Bookkeeping variables and data tracking.
    names = ["amplitude", "center", "sigma"]
    prefixes = ["g1_", "g2_", "g3_"]
    best_fit_values = [[], [], [], [], [], [], [], [], [], [], []]

    count = 0
    for i in tqdm(range(np.shape(data)[1])):
        for j in tqdm(range(np.shape(data)[2]), leave=False):
            
            # Record X and Y pixels
            best_fit_values[0].append(i)
            best_fit_values[1].append(j)
            
            # Test pixel index limiter
            if test_pixel:
                if (i, j) not in test_pixels:
                    continue 
            
            # Initialize best fit flag
            best_fit_type = "single"
            
            # Extract spaxel flux and reject routine for spaxels with data quality flags or incomplete spectra
            full_spec_1 = data[:, i, j]
            if np.isnan(full_spec_1).any():
                
                # Populate unfit values with nans
                for idx in range(2,11):
                    best_fit_values[idx].append(np.nan)
                count += 1
                #print((count / (np.shape(data)[1] * np.shape(data)[2])) * 100)
                continue
            
            # Trim spectrum and subtract background
            trim_wavelengths, trim_flux = trim_spec(wavelengths, full_spec_1, line_dict[name][3][0], line_dict[name][3][1])
            new_wave, new_flux, transition_idx = cut_line(trim_wavelengths, trim_flux, line_dict[name][2][0], line_dict[name][2][1])
            popt, pcov = fit_continuum(new_wave, new_flux)
            bkg_sub_fluxes = trim_flux - OneDPolynomial(trim_wavelengths, popt[0], popt[1])
            
            #### Fit Single Gaussian Model ####
            
            single_g1 = GaussianModel(prefix="g1_")
            single_params = single_g1.guess(bkg_sub_fluxes, x=trim_wavelengths)
            single_result = single_g1.fit(bkg_sub_fluxes, single_params, x=trim_wavelengths)
            single_redchi = single_result.redchi
            
            
            #### Fit Two Gaussian Models ####
            
            narrow_sigma = line_dict[name][5][2] / 2.3548200
            narrow_sigma2 = line_dict[name][5][3] / 2.3548200
            broad_sigma = line_dict[name][5][3] / 2.3548200
            narrow_broad_limit = (narrow_sigma + broad_sigma)/2
            line_offset = line_dict[name][5][0]
            c1_factor = line_dict[name][5][1]
            
            #### GUESS 1 (Center + RED)
            
            double_g1_1 = GaussianModel(prefix="g1_")
            double_g2_1 = GaussianModel(prefix="g2_")
            double_params_1 = double_g1_1.guess(bkg_sub_fluxes, x=trim_wavelengths)
            
            # Center
            double_params_1.update(double_g1_1.make_params(center=dict(value=line_dict[name][4], min=line_dict[name][4] - line_offset*c1_factor, max=line_dict[name][4] + line_offset*c1_factor),
                                           sigma=dict(value=narrow_sigma, min=narrow_sigma/2, max=narrow_broad_limit),
                                           amplitude=dict(value=np.max(bkg_sub_fluxes), min=0, max=100)))
            # Red
            double_params_1.update(double_g2_1.make_params(center=dict(value=line_dict[name][4] + line_offset, min=line_dict[name][4] + line_offset*c1_factor, max=line_dict[name][4] + line_offset*3),
                                           sigma=dict(value=narrow_sigma, min=narrow_sigma/2, max=narrow_broad_limit),
                                           amplitude=dict(value=np.max(bkg_sub_fluxes) * 0.8, min=0, max=50)))
            
            double_model_1 = double_g1_1 + double_g2_1
            #double_init_1 = double_model_1.eval(double_params_1, x=trim_wavelengths)
            double_result_1 = double_model_1.fit(bkg_sub_fluxes, double_params_1, x=trim_wavelengths)
            
            #### GUESS 2 (Center + BLUE)
            
            double_g1_2 = GaussianModel(prefix="g1_")
            double_g2_2 = GaussianModel(prefix="g2_")
            double_params_2 = double_g1_2.guess(bkg_sub_fluxes, x=trim_wavelengths)
            # Center
            double_params_2.update(double_g1_2.make_params(center=dict(value=line_dict[name][4], min=line_dict[name][4] - line_offset*c1_factor, max=line_dict[name][4] + line_offset*c1_factor),
                                           sigma=dict(value=narrow_sigma, min=narrow_sigma/2, max=narrow_broad_limit),
                                           amplitude=dict(value=np.max(bkg_sub_fluxes), min=0, max=100)))
            # Blue
            double_params_2.update(double_g2_2.make_params(center=dict(value=line_dict[name][4] - line_offset, max=line_dict[name][4] - line_offset*c1_factor, min=line_dict[name][4] - line_offset*3),
                                           sigma=dict(value=narrow_sigma, min=narrow_sigma/2, max=narrow_broad_limit),
                                           amplitude=dict(value=np.max(bkg_sub_fluxes) * 0.8, min=0, max=50)))
            double_model_2 = double_g1_2 + double_g2_2 
            #double_init_2 = double_model_2.eval(double_params_2, x=trim_wavelengths)
            double_result_2 = double_model_2.fit(bkg_sub_fluxes, double_params_2, x=trim_wavelengths)
            
            #### GUESS 3 (CENTER + WEAK BLUE)
            
            double_g1_3 = GaussianModel(prefix="g1_")
            double_g2_3 = GaussianModel(prefix="g2_")
            double_params_3 = double_g1_3.guess(bkg_sub_fluxes, x=trim_wavelengths)
            # Center
            double_params_3.update(double_g1_3.make_params(center=dict(value=line_dict[name][4], min=line_dict[name][4] - line_offset*c1_factor, max=line_dict[name][4] + line_offset*c1_factor),
                                           sigma=dict(value=narrow_sigma2, min=narrow_sigma2/2, max=narrow_broad_limit),
                                           amplitude=dict(value=np.max(bkg_sub_fluxes) * 0.7, min=0, max=100)))
            # Blue
            double_params_3.update(double_g2_3.make_params(center=dict(value=line_dict[name][4] - line_offset, max=line_dict[name][4] - line_offset*c1_factor, min=line_dict[name][4] - line_offset*3),
                                           sigma=dict(value=broad_sigma, min=narrow_broad_limit, max=broad_sigma*2),
                                           amplitude=dict(value=np.max(bkg_sub_fluxes) * 0.7, min=0, max=50)))
            double_model_3 = double_g1_3 + double_g2_3 
            #double_init_3 = double_model_3.eval(double_params_3, x=trim_wavelengths)
            double_result_3 = double_model_3.fit(bkg_sub_fluxes, double_params_3, x=trim_wavelengths)
            
            #### GUESS 4 (CENTER + WEAK RED)
            
            double_g1_4 = GaussianModel(prefix="g1_")
            double_g2_4 = GaussianModel(prefix="g2_")
            double_params_4 = double_g1_4.guess(bkg_sub_fluxes, x=trim_wavelengths)
            # Center
            double_params_4.update(double_g1_4.make_params(center=dict(value=line_dict[name][4], min=line_dict[name][4] - line_offset*c1_factor, max=line_dict[name][4] + line_offset*c1_factor),
                                           sigma=dict(value=narrow_sigma2, min=narrow_sigma2/2, max=narrow_broad_limit),
                                           amplitude=dict(value=np.max(bkg_sub_fluxes) * 0.7, min=0, max=100)))
            # Red
            double_params_4.update(double_g2_4.make_params(center=dict(value=line_dict[name][4] + line_offset, min=line_dict[name][4] + line_offset*c1_factor, max=line_dict[name][4] + line_offset*3),
                                           sigma=dict(value=broad_sigma, min=narrow_broad_limit, max=broad_sigma*2),
                                           amplitude=dict(value=np.max(bkg_sub_fluxes) * 0.7, min=0, max=75)))
            double_model_4 = double_g1_4 + double_g2_4 
            #double_init_4 = double_model_4.eval(double_params_4, x=trim_wavelengths)
            double_result_4 = double_model_4.fit(bkg_sub_fluxes, double_params_4, x=trim_wavelengths)
            
            #### GUESS 5 (CENTER + BROAD CENTER)
            
            double_g1_5 = GaussianModel(prefix="g1_")
            double_g2_5 = GaussianModel(prefix="g2_")
            double_params_5 = double_g1_5.guess(bkg_sub_fluxes, x=trim_wavelengths)
            # Center
            double_params_5.update(double_g1_5.make_params(center=dict(value=line_dict[name][4], min=line_dict[name][4] - line_offset*c1_factor, max=line_dict[name][4] + line_offset*c1_factor),
                                           sigma=dict(value=narrow_sigma, min=narrow_sigma/2, max=narrow_broad_limit*2),
                                           amplitude=dict(value=np.max(bkg_sub_fluxes) * 0.7, min=0, max=1000)))
            # Red
            double_params_5.update(double_g2_5.make_params(center=dict(value=line_dict[name][4] + line_offset, min=line_dict[name][4] + line_offset*c1_factor, max=line_dict[name][4] + line_offset*3),
                                           sigma=dict(value=broad_sigma*2, min=narrow_broad_limit, max=broad_sigma*3),
                                           amplitude=dict(value=np.max(bkg_sub_fluxes) * 0.7, min=0, max=500)))
            double_model_5 = double_g1_5 + double_g2_5
            #double_init_4 = double_model_4.eval(double_params_4, x=trim_wavelengths)
            double_result_5 = double_model_5.fit(bkg_sub_fluxes, double_params_5, x=trim_wavelengths)
            
            # Record reduced chi squared values
            try:
                double_redchi_1 = double_result_1.redchi
            except:
                double_redchi_1 = 1e7
            try:
                double_redchi_2 = double_result_2.redchi
            except:
                double_redchi_2 = 1e7
            try:
                double_redchi_3 = double_result_3.redchi
            except:
                double_redchi_3 = 1e7
            try:
                double_redchi_4 = double_result_4.redchi
            except:
                double_redchi_4 = 1e7
            try:
                double_redchi_5 = double_result_5.redchi
            except:
                double_redchi_5 = 1e7
            
            double_results = [double_result_1, double_result_2, double_result_3, double_result_4, double_result_5]
            double_redchis = [double_redchi_1, double_redchi_2, double_redchi_3, double_redchi_4, double_redchi_5]
            
            double_redchi = np.min(double_redchis)
            double_result = double_results[np.argmin(double_redchis)]
            
            #### Fit Three Gaussian Models ####
            
            #### GUESS 1 (CENTER + RED + BROAD)
            
            #print(test_amp)
            
            
            
            triple_model_g1_1 = GaussianModel(prefix="g1_")
            triple_model_g2_1 = GaussianModel(prefix="g2_")
            triple_model_g3_1 = GaussianModel(prefix="g3_")
            triple_params_1 = triple_model_g1_1.guess(bkg_sub_fluxes, x=trim_wavelengths)
            test_amp = triple_params_1.valuesdict()['g1_amplitude']
            #print(test_amp)
            triple_params_1.update(triple_model_g1_1.make_params(center=dict(value=line_dict[name][4], min=line_dict[name][4] - line_offset*c1_factor, max=line_dict[name][4] + line_offset*c1_factor),
                                           sigma=dict(value=narrow_sigma, min=narrow_sigma/2, max=narrow_broad_limit),
                                           amplitude=dict(value=test_amp, min=0, max=100)))
            
            triple_params_1.update(triple_model_g2_1.make_params(center=dict(value=line_dict[name][4] + line_offset, min=line_dict[name][4] + line_offset*c1_factor, max=line_dict[name][4] + line_offset*1.3),
                                           sigma=dict(value=narrow_sigma2, min=narrow_sigma2/2, max=narrow_broad_limit),
                                           amplitude=dict(value=test_amp, min=0, max=50)))
            
            triple_params_1.update(triple_model_g3_1.make_params(center=dict(value=line_dict[name][4], min=line_dict[name][4] - line_offset*0.75, max=line_dict[name][4] + line_offset*0.75),
                                           sigma=dict(value=broad_sigma, min=narrow_broad_limit, max=broad_sigma*2),
                                           amplitude=dict(value=test_amp, min=0, max=50)))
            triple_model_1 = triple_model_g1_1 + triple_model_g2_1 + triple_model_g3_1
            triple_init_1 = triple_model_1.eval(triple_params_1, x=trim_wavelengths)
            triple_result_1 = triple_model_1.fit(bkg_sub_fluxes, triple_params_1, x=trim_wavelengths)
            
            #### GUESS 2 (CENTER + BLUE + BROAD)
            
            triple_model_g1_2 = GaussianModel(prefix="g1_")
            triple_model_g2_2 = GaussianModel(prefix="g2_")
            triple_model_g3_2 = GaussianModel(prefix="g3_")
            triple_params_2 = triple_model_g1_2.guess(bkg_sub_fluxes, x=trim_wavelengths)
            triple_params_2.update(triple_model_g1_2.make_params(center=dict(value=line_dict[name][4], min=line_dict[name][4] - line_offset*c1_factor, max=line_dict[name][4] + line_offset*c1_factor),
                                           sigma=dict(value=narrow_sigma, min=narrow_sigma/2, max=narrow_broad_limit),
                                           amplitude=dict(value=test_amp, min=0, max=100)))
            
            triple_params_2.update(triple_model_g2_2.make_params(center=dict(value=line_dict[name][4] - line_offset, max=line_dict[name][4] - line_offset*c1_factor, min=line_dict[name][4] - line_offset*1.3),
                                           sigma=dict(value=narrow_sigma2, min=narrow_sigma2/2, max=narrow_broad_limit),
                                           amplitude=dict(value=test_amp, min=0, max=50)))
            
            triple_params_2.update(triple_model_g3_2.make_params(center=dict(value=line_dict[name][4], min=line_dict[name][4] - line_offset*0.75, max=line_dict[name][4] + line_offset*0.75),
                                           sigma=dict(value=broad_sigma, min=narrow_broad_limit, max=broad_sigma*2),
                                           amplitude=dict(value=test_amp, min=0, max=50)))
            triple_model_2 = triple_model_g1_2 + triple_model_g2_2 + triple_model_g3_2
            triple_init_2 = triple_model_2.eval(triple_params_2, x=trim_wavelengths)
            triple_result_2 = triple_model_2.fit(bkg_sub_fluxes, triple_params_2, x=trim_wavelengths)
            
            
            #### GUESS 3 (CENTER + WEAK BLUE + BROAD)
            
            triple_model_g1_3 = GaussianModel(prefix="g1_")
            triple_model_g2_3 = GaussianModel(prefix="g2_")
            triple_model_g3_3 = GaussianModel(prefix="g3_")
            triple_params_3 = triple_model_g1_3.guess(bkg_sub_fluxes, x=trim_wavelengths)
            triple_params_3.update(triple_model_g1_3.make_params(center=dict(value=line_dict[name][4], min=line_dict[name][4] - line_offset*c1_factor, max=line_dict[name][4] + line_offset*c1_factor),
                                           sigma=dict(value=narrow_sigma, min=narrow_sigma/2, max=narrow_broad_limit),
                                           amplitude=dict(value=test_amp, min=0, max=100)))
            
            triple_params_3.update(triple_model_g2_3.make_params(center=dict(value=line_dict[name][4] - line_offset, max=line_dict[name][4] - line_offset*c1_factor, min=line_dict[name][4] - line_offset*1.3),
                                           sigma=dict(value=broad_sigma, min=narrow_broad_limit, max=broad_sigma*2),
                                           amplitude=dict(value=test_amp, min=0, max=50)))
            
            triple_params_3.update(triple_model_g3_3.make_params(center=dict(value=line_dict[name][4], min=line_dict[name][4] - line_offset, max=line_dict[name][4] + line_offset),
                                           sigma=dict(value=broad_sigma, min=narrow_broad_limit, max=broad_sigma*2),
                                           amplitude=dict(value=test_amp, min=0, max=50)))
            triple_model_3 = triple_model_g1_3 + triple_model_g2_3 + triple_model_g3_3
            triple_init_3 = triple_model_3.eval(triple_params_3, x=trim_wavelengths)
            triple_result_3 = triple_model_3.fit(bkg_sub_fluxes, triple_params_3, x=trim_wavelengths)
            
            #### GUESS 4 (CENTER + WEAK RED + BROAD)
            
            triple_model_g1_4 = GaussianModel(prefix="g1_")
            triple_model_g2_4 = GaussianModel(prefix="g2_")
            triple_model_g3_4 = GaussianModel(prefix="g3_")
            triple_params_4 = triple_model_g1_4.guess(bkg_sub_fluxes, x=trim_wavelengths)
            triple_params_4.update(triple_model_g1_4.make_params(center=dict(value=line_dict[name][4], min=line_dict[name][4] - line_offset*c1_factor, max=line_dict[name][4] + line_offset*c1_factor),
                                           sigma=dict(value=narrow_sigma, min=narrow_sigma/2, max=narrow_broad_limit),
                                           amplitude=dict(value=test_amp, min=0, max=100)))
            
            triple_params_4.update(triple_model_g2_4.make_params(center=dict(value=line_dict[name][4] + line_offset, min=line_dict[name][4] + line_offset*c1_factor, max=line_dict[name][4] + line_offset*1.3),
                                           sigma=dict(value=broad_sigma, min=narrow_broad_limit, max=broad_sigma*2),
                                           amplitude=dict(value=test_amp, min=0, max=50)))
            
            triple_params_4.update(triple_model_g3_4.make_params(center=dict(value=line_dict[name][4], min=line_dict[name][4] - line_offset, max=line_dict[name][4] + line_offset),
                                           sigma=dict(value=broad_sigma, min=narrow_broad_limit, max=broad_sigma*2),
                                           amplitude=dict(value=test_amp, min=0, max=50)))
            triple_model_4 = triple_model_g1_4 + triple_model_g2_4 + triple_model_g3_4
            triple_init_4 = triple_model_4.eval(triple_params_4, x=trim_wavelengths)
            triple_result_4 = triple_model_4.fit(bkg_sub_fluxes, triple_params_4, x=trim_wavelengths)
            
            #### GUESS 5 (CENTER + Double BROAD)
            
            triple_model_g1_5 = GaussianModel(prefix="g1_")
            triple_model_g2_5 = GaussianModel(prefix="g2_")
            triple_model_g3_5 = GaussianModel(prefix="g3_")
            triple_params_5 = triple_model_g1_5.guess(bkg_sub_fluxes, x=trim_wavelengths)
            triple_params_5.update(triple_model_g1_5.make_params(center=dict(value=line_dict[name][4], min=line_dict[name][4] - line_offset*c1_factor, max=line_dict[name][4] + line_offset*c1_factor),
                                           sigma=dict(value=narrow_sigma, min=narrow_sigma/2, max=narrow_broad_limit*2),
                                           amplitude=dict(value=test_amp, min=0, max=1000)))
            
            triple_params_5.update(triple_model_g2_5.make_params(center=dict(value=line_dict[name][4], min=line_dict[name][4] - line_offset*c1_factor, max=line_dict[name][4] + line_offset*c1_factor),
                                           sigma=dict(value=narrow_sigma, min=narrow_broad_limit, max=broad_sigma*3),
                                           amplitude=dict(value=test_amp, min=0, max=500)))
            
            triple_params_5.update(triple_model_g3_5.make_params(center=dict(value=line_dict[name][4], min=line_dict[name][4] - line_offset*c1_factor, max=line_dict[name][4] + line_offset*c1_factor),
                                           sigma=dict(value=broad_sigma*2, min=narrow_broad_limit, max=broad_sigma*3),
                                           amplitude=dict(value=test_amp, min=0, max=400)))
            triple_model_5 = triple_model_g1_5 + triple_model_g2_5 + triple_model_g3_5
            #triple_init_5 = triple_model_5.eval(triple_params_5, x=trim_wavelengths)
            triple_result_5 = triple_model_5.fit(bkg_sub_fluxes, triple_params_5, x=trim_wavelengths)
            
            
            
            # Record reduced chi squared values
            try:
                triple_redchi_1 = triple_result_1.redchi
            except:
                triple_redchi_1 = 1e7
            try:
                triple_redchi_2 = triple_result_2.redchi
            except:
                triple_redchi_2 = 1e7
            try:
                triple_redchi_3 = triple_result_3.redchi
            except:
                triple_redchi_3 = 1e7
            try:
                triple_redchi_4 = triple_result_4.redchi
            except:
                triple_redchi_4 = 1e7
            try:
                triple_redchi_5 = triple_result_5.redchi
            except:
                triple_redchi_5 = 1e7
            
            triple_results = [triple_result_1, triple_result_2, triple_result_3, triple_result_4, triple_result_5]
            triple_redchis = [triple_redchi_1, triple_redchi_2, triple_redchi_3, triple_redchi_4, triple_redchi_5]
            
            triple_redchi = np.min(triple_redchis)
            triple_result = triple_results[np.argmin(triple_redchis)]
        
            
            # Evaluate which model is better
            
            
            if single_redchi < double_redchi:
                if single_redchi < triple_redchi:
                    base_array[i][j] = 1
                    lowest_redchi = single_redchi
                    base_array2[i][j] = lowest_redchi
                    best_fit_type = "single"
                else:
                    base_array[i][j] = 3
                    lowest_redchi = triple_redchi
                    base_array2[i][j] = lowest_redchi
                    best_fit_type = "triple"
            else:
                if double_redchi < triple_redchi:
                    base_array[i][j] = 2
                    lowest_redchi = double_redchi
                    base_array2[i][j] = lowest_redchi
                    best_fit_type = "double"
                else:
                    base_array[i][j] = 3
                    lowest_redchi = triple_redchi
                    base_array2[i][j] = lowest_redchi
                    best_fit_type = "triple"
            
            if best_fit_type == "single":
                best_result = single_result                
            if best_fit_type == "double":
                best_result = double_result
            if best_fit_type == "triple":
                best_result = triple_result

            # Record best values
            best_idx = 2
            for prefix in prefixes:
                for pname in names:
                    if prefix+pname in best_result.best_values:
                        best_fit_values[best_idx].append(best_result.best_values[prefix+pname])
                        best_idx += 1
                    else:
                        best_fit_values[best_idx].append(np.nan)
                        best_idx += 1
            
            
            if test_pixel:
                
                fig, ax = plt.subplots()
                peaks, _ = find_peaks(bkg_sub_fluxes, height=np.max(bkg_sub_fluxes)*0.2)
                trim_peaks = []
                for peak in peaks:
                    if trim_wavelengths[peak] > line_dict[name][2][0] and trim_wavelengths[peak] < line_dict[name][2][1]:
                        trim_peaks.append(peak)
                ax.scatter(trim_wavelengths[trim_peaks], bkg_sub_fluxes[trim_peaks], marker="*", s=100, c="gold", label=f"PEAKS: n={len(trim_peaks)}", zorder=0)
                ax.scatter(trim_wavelengths, bkg_sub_fluxes, c="white", alpha=0.75)
                ax.plot(trim_wavelengths, single_result.best_fit, c="cyan", label=f"One Gaussian: {round(single_redchi, 2)}", lw=0.9)
                ax.plot(trim_wavelengths, double_result.best_fit, c="lime", label=f"Two Gaussian: {round(double_redchi, 2)}", lw=0.9)
                ax.plot(trim_wavelengths, triple_result.best_fit, c="orangered", label=f"Three Gaussian: {round(triple_redchi, 2)}", lw=0.9)
                #print(double_result.fit_report())
                #print(triple_result.fit_report())
                ax.legend()
                plt.savefig(f"./../diagnostic_plots/dynamic_multicomponent/{name}/test_spaxels/{j}_{i}.pdf", bbox_inches="tight")
                plt.close()
            count += 1
            #print((count / (np.shape(data)[1] * np.shape(data)[2])) * 100)

    fitparams = Table(best_fit_values, names=("XPIX", "YPIX", "G1AMP", "G1CEN", "G1SIGMA", "G2AMP", "G2CEN", "G2SIGMA", "G3AMP", "G3CEN", "G3SIGMA"))
    fitparams.write(f"./../diagnostic_plots/dynamic_multicomponent/{name}/fit.dat", format="ipac", overwrite=True)

core = "S"
name = "[FeII]"
type = "triple"
multicomponent = True
redshift = 0.044601
main_routine(core, name, redshift, test_pixel=False, experimental_bkg=False)#, test_index=2)







