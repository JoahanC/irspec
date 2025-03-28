
import numpy as np 
import astropy.units as u
from irspec.spec_helpers import trim_spec, cut_line, fit_continuum, find_nans
from astropy.io import ascii
from specutils import Spectrum1D, SpectralRegion
from specutils.manipulation import extract_region
from specutils.analysis import snr_derived

from irspec.fitfuncs import * 
from lmfit.models import GaussianModel
from lmfit import Parameters
from astropy.table import Table
from tqdm import tqdm

import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.constants as const
from matplotlib.colors import LogNorm
from astropy.visualization.wcsaxes import add_scalebar
import matplotlib.font_manager as fm

from irspec.plotparams import PlotParams
from irspec.emission_io import read_line_params
pltparams = PlotParams(scaling="presentation")

MODES = ["cube", "spaxel", "spaxels"]
TEST_SPAXELS = [(32, 10), (26, 16), (31, 19), (26, 36), (25, 21), (28, 38), 
                (25, 22), (35, 12), (29, 10), (29, 32), (22, 26), (26, 22),
                (19, 26), (26, 19), (31, 14)]

class SpaxelFit:
    """ 
    A class for performing multicomponent gaussian fits to IFS data
    """
    
    def __init__(self, datacube, name, output, 
                 mode="cube", test_spaxel=(20,20), test_spaxels=None):
        """ 
        Arguments
        ---------
        datacube : Datacube
            A Datacube metadata container.
        name : str
            The name of the emission line that is to be fit.
        redshift : float
            The redshift corresponding to the datacube.
        output : str
            The output file name to write out the fitting results 
        mode : str
            This should take on either `cube`, `spaxel` or `spaxels`. `cube` is the
            default and fits the entire datacube. `spaxel` fits a designated 
            spaxel. `spaxels` fits a list of predetermined spaxels.
        test_spaxel : tuple
            This is the coordinates of the spaxel to fit in `spaxel` mode. The 
            default is (20, 20).
        test_spaxels : list
            A list of tuples corresponding to the spaxels to trial in `spaxels` mode
        """
        
        if mode not in MODES:
            raise ValueError("Mode must be 'cube', 'spaxel', or 'spaxels'")
        if test_spaxels == None:
            self.test_spaxels = TEST_SPAXELS
        
        self.datacube = datacube 
        self.name = name 
        self.output = output
        self.mode = mode
        self.test_spaxel = test_spaxel
        self.line_dict = read_line_params()       
        self.label = None
    
    
    def load_fit(self, filepath):
        """Loads a saved instance of SpaxelFit"""
        self.fitparams = ascii.read(filepath, format="ipac")  
    
    
    def _mask_spaxel(self, spaxel):
        """Develops the masks necessary to perform gaussian fitting."""
        #lower_continuum_region = SpectralRegion(np.min(self.datacube.wvs.value) * self.datacube._wv_unit,
        #                                        self.line_dict[self.name][2][0] * self.datacube._wv_unit)
        #upper_continuum_region = SpectralRegion(self.line_dict[self.name][2][1] * self.datacube._wv_unit,
        #                                        np.max(self.datacube.wvs.value) * self.datacube._wv_unit)
        lower_continuum_region = SpectralRegion(self.line_dict[self.name][2][0] * self.datacube._wv_unit,
                                                self.line_dict[self.name][3][0] * self.datacube._wv_unit)
        upper_continuum_region = SpectralRegion(self.line_dict[self.name][2][1] * self.datacube._wv_unit,
                                                self.line_dict[self.name][3][1] * self.datacube._wv_unit)
        continuum_region = lower_continuum_region + upper_continuum_region
        line_region = SpectralRegion(self.line_dict[self.name][3][0] * self.datacube._wv_unit,
                                        self.line_dict[self.name][3][1] * self.datacube._wv_unit)
        spaxel_continuum = extract_region(spaxel, continuum_region, return_single_spectrum=True)
        spaxel_line = extract_region(spaxel, line_region)
        return spaxel_continuum, spaxel_line


    def background_subtraction(self, spaxel_continuum, spaxel_line):
        """Performs background subtraction on a spaxel."""
        popt, pcov = fit_continuum(spaxel_continuum.spectral_axis.value, 
                                   spaxel_continuum.flux.value)
        spaxel_continuum_fit = Spectrum1D(OneDPolynomial(spaxel_line.spectral_axis.value, popt[0], popt[1]) * u.Jy, 
                                            spaxel_line.spectral_axis)
        spaxel_continuum_sub = spaxel_line - spaxel_continuum_fit
        line_fluxes = spaxel_continuum_sub.flux.value
        line_wavelengths = spaxel_continuum_sub.spectral_axis.value
        return line_wavelengths, line_fluxes, spaxel_continuum_fit.flux.value
    
    
    def two_gaussian_fit(self):
        """ 
        A simple fitting routine which fits two gaussian 
        components to a line.
        """
        
        # Bookkeeping variables and data tracking.
        names = ["amplitude", "center", "sigma"]
        prefixes = ["g1_", "g2_", "g3_"]
        best_fit_values = [[], [], [], [], [], [], [], [], [], [], [], [], []]
        
        # Iterate through all spaxels
        for y_pix in tqdm(range(self.datacube.im_shape[0])):
            for x_pix in tqdm(range(self.datacube.im_shape[0]), leave=False):
                
                # Record X and Y pixels
                best_fit_values[0].append(y_pix)
                best_fit_values[1].append(x_pix)
                
                # Test for mode case
                if self.mode == "spaxels":
                    if (x_pix, y_pix) not in self.test_spaxels:
                        continue
                if self.mode == "spaxel":
                    if (x_pix, y_pix) != self.test_spaxel:
                        continue
                
                (flux, flux_err, dq) = self.datacube.spaxel_values(x_pix, y_pix)
                if 513 in np.unique(dq):
                    for idx in range(2,13):
                        best_fit_values[idx].append(np.nan)
                    continue
                spaxel = Spectrum1D(flux, self.datacube.wvs)
                
                # Process spaxel
                spaxel_continuum, spaxel_line = self._mask_spaxel(spaxel)
                line_wavelengths, line_fluxes, _ = self.background_subtraction(spaxel_continuum, spaxel_line)
                
                # Estimate line parameters
                line_center = self.line_dict[self.name][4]
                center_offset = self.line_dict[self.name][5]
                center_cutoff = self.line_dict[self.name][6]
                narrow_sigma = self.line_dict[self.name][7] / 2.3548200
                broad_sigma = self.line_dict[self.name][8] / 2.3548200
                narrow_amp = np.sqrt(2*np.pi)*np.max(line_fluxes)*narrow_sigma
                broad_amp = np.sqrt(2*np.pi)*np.max(line_fluxes)*broad_sigma
                
                ### Run lmfit fitting routine ###
                
                # Narrow center + narrow red
                double_g1_2 = GaussianModel(prefix="g1_")
                double_g2_2 = GaussianModel(prefix="g2_")
                double_params_2 = double_g1_2.guess(line_fluxes, x=line_wavelengths)
                double_params_2.update(double_g1_2.make_params(center=dict(value=line_center, min=line_center-center_offset, max=line_center+center_offset),
                                            sigma=dict(value=narrow_sigma),
                                            amplitude=dict(value=narrow_amp, max=narrow_amp*3, min=narrow_amp/2)))
                double_params_2.update(double_g2_2.make_params(center=dict(value=line_center+center_offset, min=line_center, max=line_center+center_cutoff),
                                            sigma=dict(value=narrow_sigma),
                                            amplitude=dict(value=narrow_amp, max=narrow_amp*2, min=narrow_amp/2)))
                double_model_2 = double_g1_2 + double_g2_2
                double_result_2 = double_model_2.fit(line_fluxes, double_params_2, x=line_wavelengths)
                
                # Narrow center + narrow blue
                double_g1_3 = GaussianModel(prefix="g1_")
                double_g2_3 = GaussianModel(prefix="g2_")
                double_params_3 = double_g1_3.guess(line_fluxes, x=line_wavelengths)
                double_params_3.update(double_g1_3.make_params(center=dict(value=line_center, min=line_center-center_offset, max=line_center+center_offset),
                                            sigma=dict(value=narrow_sigma),
                                            amplitude=dict(value=narrow_amp, max=narrow_amp*3, min=narrow_amp/2)))
                double_params_3.update(double_g2_3.make_params(center=dict(value=line_center-center_offset, min=line_center-center_cutoff, max=line_center),
                                            sigma=dict(value=narrow_sigma),
                                            amplitude=dict(value=narrow_amp, max=narrow_amp*2, min=narrow_amp/2)))
                double_model_3 = double_g1_3 + double_g2_3
                double_result_3 = double_model_3.fit(line_fluxes, double_params_3, x=line_wavelengths)
                
                # Narrow center + broad red
                double_g1_4 = GaussianModel(prefix="g1_")
                double_g2_4 = GaussianModel(prefix="g2_")
                double_params_4 = double_g1_4.guess(line_fluxes, x=line_wavelengths)
                double_params_4.update(double_g1_4.make_params(center=dict(value=line_center, min=line_center-center_offset, max=line_center+center_offset),
                                            sigma=dict(value=narrow_sigma),
                                            amplitude=dict(value=narrow_amp, max=narrow_amp*3, min=narrow_amp/2)))
                double_params_4.update(double_g2_4.make_params(center=dict(value=line_center+center_offset, min=line_center, max=line_center+center_cutoff),
                                            sigma=dict(value=broad_sigma),
                                            amplitude=dict(value=broad_amp, max=broad_amp*2, min=broad_amp/2)))
                double_model_4 = double_g1_4 + double_g2_4
                double_result_4 = double_model_4.fit(line_fluxes, double_params_4, x=line_wavelengths)
                
                # Narrow center + broad blue
                double_g1_5 = GaussianModel(prefix="g1_")
                double_g2_5 = GaussianModel(prefix="g2_")
                double_params_5 = double_g1_5.guess(line_fluxes, x=line_wavelengths)
                double_params_5.update(double_g1_5.make_params(center=dict(value=line_center, min=line_center-center_offset, max=line_center+center_offset),
                                            sigma=dict(value=narrow_sigma),
                                            amplitude=dict(value=narrow_amp, max=narrow_amp*3, min=narrow_amp/2)))
                double_params_5.update(double_g2_5.make_params(center=dict(value=line_center-center_offset, min=line_center-center_cutoff, max=line_center),
                                            sigma=dict(value=broad_sigma),
                                            amplitude=dict(value=broad_amp, max=broad_amp*2, min=broad_amp/2)))
                double_model_5 = double_g1_5 + double_g2_5
                double_result_5 = double_model_5.fit(line_fluxes, double_params_2, x=line_wavelengths)
                
                # Record reduced chi squared values
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
                
                double_results = [double_result_2, double_result_3, double_result_4, double_result_5]
                double_redchis = [double_redchi_2, double_redchi_3, double_redchi_4, double_redchi_5]
                
                double_redchi = np.min(double_redchis)
                double_result = double_results[np.argmin(double_redchis)]
            
                best_idx = 2
                for prefix in prefixes:
                    for pname in names:
                        if prefix+pname in double_result.best_values:
                            best_fit_values[best_idx].append(double_result.best_values[prefix+pname])
                            best_idx += 1
                        else:
                            best_fit_values[best_idx].append(np.nan)
                            best_idx += 1
                
                best_fit_values[12].append(2)
                best_fit_values[11].append(double_redchi)
                
                
                
                
        if self.mode == "cube":
            for listy in best_fit_values:
                print(len(listy))
            self.fitparams = Table(best_fit_values, names=("XPIX", "YPIX", "G1AMP", "G1CEN", "G1SIGMA", "G2AMP", "G2CEN", "G2SIGMA", "G3AMP", "G3CEN", "G3SIGMA", "REDCHI", "NCOMP"))
            self.fitparams.write(self.output + "twogaussian_raw.dat", format="ipac", overwrite=True)
            self.label = "twogaussian"
    
    
    def render_spaxel_fit(self, x_pix, y_pix):
        """Renders the multicomponent gaussian fit for an individual spaxel."""
        (flux, flux_err, dq) = self.datacube.spaxel_values(x_pix, y_pix)
        spaxel = Spectrum1D(flux, self.datacube.wvs)
        spaxel_continuum, spaxel_line = self._mask_spaxel(spaxel)
        line_wavelengths, line_fluxes, background_flux = self.background_subtraction(spaxel_continuum, spaxel_line)
        # Generate component fluxes and total flux
        model_components = []
        total_flux = np.copy(background_flux)
        spaxel_idx = x_pix * self.datacube.im_shape[0] + y_pix
        for idx in range(1, int(self.fitparams["NCOMP"][spaxel_idx]) + 1):
            params = Parameters()
            params.add(name="amplitude", value=self.fitparams[f"G{idx}AMP"][spaxel_idx])
            params.add(name="center", value=self.fitparams[f"G{idx}CEN"][spaxel_idx])
            params.add(name="sigma", value=self.fitparams[f"G{idx}SIGMA"][spaxel_idx])
            model_component = GaussianModel()
            component_flux = model_component.eval(params=params, x=line_wavelengths)
            total_flux += component_flux
            model_components.append(background_flux + component_flux)
        
        # Plot
        fig, ax = plt.subplots()
        fig.set_size_inches(10, 10)
        ax.scatter(spaxel_line.spectral_axis.value, spaxel_line.flux.value, label="Data", c="white")
        ax.plot(spaxel_line.spectral_axis.value, background_flux, label="Continuum", c="yellow", ls="dashed")
        ax.plot(spaxel_line.spectral_axis.value, total_flux, label="Model", c="cyan")
        for idx, component in enumerate(model_components):
            ax.plot(spaxel_line.spectral_axis.value, component, label=f"Component {idx+1}")
        ax.set_title(f"Spaxel: ({x_pix}, {y_pix})", loc="right")
        ax.set_title(f"Gaussian Decomposition", loc="left")
        ax.set_xlabel("Wavelength [nm]")
        ax.set_ylabel("Flux [Jy]")
        ax.legend()
        plt.show()
        
        
        
    """def render_current_fit(self):
        
        keys = list(self.fitparams.columns)[2:]
        pprams = PlotParams(palatte="dark", scaling="presentation")
        render_vel_disp_rel_vel_scatter(self.fitparams, self.line_dict[self.name][4], self.name, self.output, savefig=True)
        for key in keys:
            if "3" not in key:
                if "AMP" in key:
                    render_amplitude_plot(self.fitparams, self.name, self.line_dict[self.name][4], wcs, self.output, param=key, savefig=True)
                if "CEN" in key:
                    render_rel_vel_plot(self.fitparams, self.line_dict[self.name][4], self.name, wcs, self.output, param=key, savefig=True)
                if "SIGMA" in key:
                    render_vel_disp_plot(self.fitparams, self.line_dict[self.name][4], self.name, wcs, self.output, param=key, savefig=True)"""
    
    
    def render_multicomponent_plot(self, savefig=False):
        """Renders a spaxel map illustrating the number of components fit to 
        each spaxel."""
        
        base_array = np.zeros((np.max(self.fitparams["XPIX"]) + 1, np.max(self.fitparams["YPIX"]) + 1))
        for idx, _ in enumerate(self.fitparams["XPIX"]):
            if np.isnan(self.fitparams["G3AMP"][idx]):
                if np.isnan(self.fitparams["G2AMP"][idx]):
                    if np.isnan(self.fitparams["G1AMP"][idx]):
                        base_array[self.fitparams["XPIX"][idx]][self.fitparams["YPIX"][idx]] = 0
                    else:
                        base_array[self.fitparams["XPIX"][idx]][self.fitparams["YPIX"][idx]] = 1
                else:
                    base_array[self.fitparams["XPIX"][idx]][self.fitparams["YPIX"][idx]] = 2
            else:
                base_array[self.fitparamsa["XPIX"][idx]][self.fitparams["YPIX"][idx]] = 3
        
        fig = plt.figure()
        ax = plt.subplot()
        fig.set_size_inches(10, 8)
        cmap = plt.get_cmap('bone', np.max(base_array) - np.min(base_array) + 1)
        image = ax.imshow(base_array, cmap=cmap, vmin=np.min(base_array) - 0.5, vmax=np.max(base_array) + 0.5)
        if self.line_dict[self.name][0] == 1:
            ax.scatter(24, 28, c="white", edgecolors="black", marker="*", s=1000)
        if self.line_dict[self.name][0] == 2:
            ax.scatter(21, 26, c="white", edgecolors="black", marker="*", s=1000)
        if self.line_dict[self.name][0] == 3:
            ax.scatter(21, 25, c="white", edgecolors="black", marker="*", s=1000)
        if self.line_dict[self.name][0] == 4:
            ax.scatter(18, 21, c="white", edgecolors="black", marker="*", s=1000)
        plt.colorbar(image, ticks=np.arange(np.min(base_array), np.max(base_array) + 1))
        ax.set_xlabel("Right Ascension")
        ax.set_ylabel("Declination")
        ax.set_title("# of Components", loc="right")
        title_name = self.name
        if "H_2_S_1" in self.name:
            title_name = r"H$_{2}$ S(1)"
        if "H_2_S_3" in self.name:
            title_name = r"H$_{2}$ S(3)"
        if "14" in self.name:
            title_name = "[NeV]"+r"$_{14}$"
        ax.set_title(title_name, loc="left")
        #gc_distance = 194.99 * u.Mpc
        #scalebar_length = 5 * u.kpc
        #scalebar_angle = (scalebar_length / gc_distance).to(
        #    u.deg, equivalencies=u.dimensionless_angles()
        #)
        fontprops = fm.FontProperties(size=24, family='Helvetica')
        #add_scalebar(ax, scalebar_angle, label="5 kpc", color="white", fontproperties=fontprops)
        if savefig:
            plt.savefig(f"./../diagnostic_plots/dynamic_multicomponent/{self.name}/compnum_white.png", dpi=600)
            plt.close()
        else:
            plt.show()
    
    def render_totflux_plot(self, savefig=False):
        """
        Renders a spaxel map illustrating the number of components fit to 
        each spaxel
        """
        
        # Initialize base aarray
        base_array = np.zeros((np.max(self.fitparams["XPIX"]) + 1, np.max(self.fitparams["YPIX"]) + 1))
        
        for idx, _ in enumerate(self.fitparams["XPIX"]):
            flux_val = 0
            if not np.isnan(self.fitparams["G1AMP"][idx]):
                
                g1_gamma = np.abs(self.fitparams["G1SIGMA"][idx]) * 2.355 / self.fitparams["G1CEN"][idx]
                g1_flux = (1 / np.sqrt(np.pi * np.log(2)) * (np.pi * const.c.to('micron/s') / 2) * (self.fitparams["G1AMP"][idx] * u.Jy * g1_gamma / (self.fitparams["G1CEN"][idx] * u.micron))).to(u.watt / u.meter ** 2)
                flux_val += g1_flux.value
                
                if not np.isnan(self.fitparams["G2AMP"][idx]):
                    g2_gamma = np.abs(self.fitparams["G2SIGMA"][idx]) * 2.355 / self.fitparams["G2CEN"][idx]
                    g2_flux = (1 / np.sqrt(np.pi * np.log(2)) * (np.pi * const.c.to('micron/s') / 2) * (self.fitparams["G2AMP"][idx] * u.Jy * g2_gamma / (self.fitparams["G2CEN"][idx] * u.micron))).to(u.watt / u.meter ** 2)
                    flux_val += g2_flux.value
                    if not np.isnan(self.fitparams["G3AMP"][idx]):
                        g3_gamma = np.abs(self.fitparams["G3SIGMA"][idx]) * 2.355 / self.fitparams["G3CEN"][idx]
                        g3_flux = (1 / np.sqrt(np.pi * np.log(2)) * (np.pi * const.c.to('micron/s') / 2) * (self.fitparams["G3AMP"][idx] * u.Jy * g3_gamma / (self.fitparams["G3CEN"][idx] * u.micron))).to(u.watt / u.meter ** 2)
                        flux_val += g3_flux.value
                        if flux_val < 0:
                            base_array[self.fitparams["XPIX"][idx]][self.fitparams["YPIX"][idx]] = 0
                        else:
                            base_array[self.fitparams["XPIX"][idx]][self.fitparams["YPIX"][idx]] = flux_val
                    else:

                        if flux_val < 0:
                            base_array[self.fitparams["XPIX"][idx]][self.fitparams["YPIX"][idx]] = 0
                        else:
                            base_array[self.fitparams["XPIX"][idx]][self.fitparams["YPIX"][idx]] = flux_val
                else:

                    if flux_val < 0:
                        base_array[self.fitparams["XPIX"][idx]][self.fitparams["YPIX"][idx]] = 0
                    else:
                        base_array[self.fitparams["XPIX"][idx]][self.fitparams["YPIX"][idx]] = flux_val
            else:
                base_array[self.fitparams["XPIX"][idx]][self.fitparams["YPIX"][idx]] = flux_val

        fig = plt.figure()
        ax = plt.subplot()
        fig.set_size_inches(10, 8)
        cmap = plt.get_cmap('plasma')
        image = ax.imshow(base_array, cmap=cmap, norm=LogNorm(), origin="lower")
        if self.line_dict[self.name][0] == 1:
            ax.scatter(24, 28, c="white", edgecolors="black", marker="*", s=1000)
        if self.line_dict[self.name][0] == 2:
            ax.scatter(21, 26, c="white", edgecolors="black", marker="*", s=1000)
        if self.line_dict[self.name][0] == 3:
            ax.scatter(21, 25, c="white", edgecolors="black", marker="*", s=1000)
        if self.line_dict[self.name][0] == 4:
            ax.scatter(18, 21, c="white", edgecolors="black", marker="*", s=1000)
        cax = plt.colorbar(image)
        cax.set_label(r"[W/m$^2$]", fontsize=24, rotation=270, labelpad=25)
        ax.set_xlabel("Right Ascension")
        ax.set_ylabel("Declination")
        ax.set_title("Total Flux", loc="right")
        title_name = self.name
        if "H_2_S_1" in self.name:
            title_name = r"H$_{2}$ S(1)"
        if "H_2_S_3" in self.name:
            title_name = r"H$_{2}$ S(3)"
        if "14" in self.name:
            title_name = "[NeV]"+r"$_{14}$"
        ax.set_title(title_name, loc="left")
        #c_distance = 194.99 * u.Mpc
        #calebar_length = 5 * u.kpc
        #scalebar_angle = (scalebar_length / gc_distance).to(
        #    u.deg, equivalencies=u.dimensionless_angles()
        #)
        fontprops = fm.FontProperties(size=24, family='Helvetica')
        #add_scalebar(ax, scalebar_angle, label="5 kpc", color="white", fontproperties=fontprops)
        if savefig:
            plt.savefig(f"./../diagnostic_plots/dynamic_multicomponent/{self.name}/flux_white.png", dpi=600)
            plt.close()
        else:
            plt.show()
    
    
    def render_rel_vel_plot(self, param="G2CEN", savefig=False):
        """
        Renders a spaxel map illustrating the number of components fit to 
        each spaxel
        """
        
        # Initialize base aarray
        max_val = 0
        min_val = 0
        base_array = np.zeros((np.max(self.fitparams["XPIX"]) + 1, np.max(self.fitparams["YPIX"]) + 1))
        #all_relvels = []
        for idx, _ in enumerate(self.fitparams["XPIX"]):
            rel_vel = self.datacube.wv_to_vel(self.fitparams[param][idx], self.line_dict[self.name][4])
            base_array[self.fitparams["XPIX"][idx]][self.fitparams["YPIX"][idx]] = rel_vel.value
            #all_relvels.append(rel_vel)

        
        #high_percentile = np.nanpercentile(all_relvels, 97.5)
        
        if "1" in param:
            comp_name = r"$v_{1}$"
        if "2" in param:
            comp_name = r"$v_{2}$"
        if "3" in param:
            comp_name = r"$v_{3}$"
        
        fig = plt.figure()
        ax = plt.subplot()
        fig.set_size_inches(10, 8)
        cmap = plt.get_cmap('RdBu_r')
        image = ax.imshow(base_array, cmap=cmap, origin="lower")
        ax.scatter(21, 24, c="white", edgecolors="black", marker="*", s=1000)
        cax = plt.colorbar(image)
        cax.set_label("[km/s]", fontsize=24, rotation=270, labelpad=25)
        ax.set_xlabel("Right Ascension")
        ax.set_ylabel("Declination")
        ax.set_title(comp_name, loc="right")
        ax.set_title(self.name, loc="left")
        #gc_distance = 194.99 * u.Mpc
        #scalebar_length = 5 * u.kpc
        #scalebar_angle = (scalebar_length / gc_distance).to(
        #    u.deg, equivalencies=u.dimensionless_angles()
        #)
        fontprops = fm.FontProperties(size=24, family='Helvetica')
        #add_scalebar(ax, scalebar_angle, label="5 kpc", color="white", fontproperties=fontprops)
        if savefig:
            plt.savefig(f"./../diagnostic_plots/dynamic_multicomponent/{self.name}/{param}.png", dpi=600)
            plt.close()
        else:
            plt.show()
    
    def gaussian_integral(amplitude, center, sigma, x_1, x_2):
        return amplitude * (sp.erf((center - x_1)) / (np.sqrt(2) * sigma) - sp.erf((center - x_2)) / (np.sqrt(2) * sigma)) / 2
    
    def snr_cut(self, snr_threshold):
        """ 
        Removes and flattens fitted component array based on an SNR threshold.
        """
        
        for i in range(1, 4):
            pass
        
        pass
    
    
    def calc_line_flux(self):
        pass
    
    #def remove_component(self):
    #    
    
    