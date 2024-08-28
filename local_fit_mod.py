""" 
Helper file for CubeSpec
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import peak_widths
import astropy.units as u
import astropy.constants as const


from fitfuncs import *
from plotparams import PlotParams

from CAFE.component_model import pah_drude, gauss_prof, drude_prof, drude_int_fluxes


class LocalFit:
    
    
    def __init__(self, wave, flux, wave_range, wave_c, name, x_err=None, y_err=None):
        """ 
        A container for hosting cutout fit parameters used in localized
        fitting routines for CubeSpec.
        
        Arguments
        ---------
        wave : ndarray
            An array containing all of the wavelengths
        flux : ndarray
            An array containing all of the fluxes
        wave_range : list
            A two element list defining the wavelength cutoffs for the 
            localized fit
        wave_c : float
            The wavelength corresponding to the line to be fit.
        name : str
            The identifier name to be affixed to this fit.
        x_err : ndarray
            Any error values associated with wavelengths
        y_err : ndarray
            Any error values associated with fluxes
        """
        
        # Import and validate data
        self.wave_c = wave_c
        self.name = name
        if type(wave_range) != list:
            raise TypeError("wave_range must be a list!")
        if len(wave_range) != 2:
            raise IndexError("wave_range must have only two elements!")
        self.wave_range = wave_range
        
        #Cutout relevant portion
        min_dif_array = np.absolute(wave-wave_range[0])
        max_dif_array = np.absolute(wave-wave_range[1])
        low_idx = min_dif_array.argmin()
        high_idx = max_dif_array.argmin() 
        self.wave = wave[low_idx:high_idx]
        self.flux = flux[low_idx:high_idx]
        self.peak_idx = (np.abs(self.wave - self.wave_c)).argmin()
        
        # Define errors if provided
        if x_err == None:
            self.wave_err = np.zeros(np.shape(self.wave))
        if y_err == None:
            self.flux_err = np.zeros(np.shape(self.flux))
        
        # Define fitting fields
        self.fitflag = None
        self.fit_function = None
        self.continuum = None 
        self.gaussians = []
        
        # Define plotting style
        pltparams = PlotParams(scaling="presentation")
        self.plotcolors = pltparams.dark_colors()
        self.linestyles = pltparams.linestyles()
        self.markers = pltparams.markers()
    
    
    def main_fit(self, npoly=1, ngauss=2, spaxel_fit=False):
        """ 
        The main refitting routine. The continuum function is defined
        as an `n` order polynomial. The line is assumed to be 
        comprised of strictly gaussian components.
        
        Arguments
        ---------
        npoly : int
            The degree of the polynomial to be fit as the continuum 
            function. The default value is 1.
        ngauss : int 
            The number of gaussian components to be fit as the line 
            feature. The default value is 2.
        """
        
        # Define the initial parameter map
        self.n_params = (npoly + 1) + ngauss * 3
        self.params = [None] * self.n_params
        self.param_map = {}
        self.param_map_r = {}
        self.guess = [0] * (npoly + ngauss * 3 + 1)
        self.lower_bounds = [-np.inf] * (npoly + ngauss * 3 + 1)
        self.upper_bounds = [np.inf] * (npoly + ngauss * 3 + 1)
        self.npoly = npoly 
        self.ngauss = ngauss
        
        # Construct the parameter map
        gauss_idx = 1
        ngauss_idx = 1
        for idx, _ in enumerate(self.params):
            if idx < (npoly + 1):
                self.param_map[f"poly_{npoly - idx}"] = idx
                self.param_map_r[str(idx)] = f"poly_{npoly - idx}"
            if idx >= (npoly + 1):
                self.param_map[f"gauss_{ngauss_idx}_{gauss_idx}"] = idx
                self.param_map_r[str(idx)] = f"gauss_{ngauss_idx}_{gauss_idx}"
                gauss_idx += 1
                if gauss_idx == 4:
                    gauss_idx = 1
                    ngauss_idx += 1
        
        # Define the continuum function
        if npoly == 1:
            self.continuum = OneDPolynomial
        if npoly == 2:
            self.continuum = TwoDPolynomial
        
        # Define the complete set of gaussians
        for _ in range(ngauss):
            self.gaussians.append(OneDGaussian)
        
        
        def fitting_function(wave, *pars):
            """ 
            Dynamically defined fitting function.
            """
            if self.npoly == 1:
                continuum = self.continuum(wave, pars[0], pars[1])
            if self.npoly == 2:
                continuum = self.continuum(wave, pars[0], pars[1], pars[2])
            gaussian_values = []
            for idx, gaussian in enumerate(self.gaussians):
                gaussian_values.append(gaussian(wave, pars[npoly + 3 * idx + 1], pars[npoly + 3 * idx + 2], pars[npoly + 3 * idx + 3]))
            fluxes = continuum
            for arrays in gaussian_values:
                fluxes += arrays
            return fluxes
        
        if npoly == 1:
            self.guess[0] = (self.flux[-1] - self.flux[0]) / (self.wave[-1] - self.wave[0])
            self.guess[1] = self.flux[0] - self.guess[0] * self.wave[0]
        if npoly == 2:
            self.guess[0] = 0
            self.guess[1] = 0
            self.guess[2] = 0
        
        amp_factors = [1, 1, 0.2, 0.1]
        broad_factors = [1, 1, 2, 3]
        width_data = peak_widths(self.flux, [self.peak_idx])

        for idx in range(ngauss):
            self.guess[npoly + 3 * idx + 1] = np.max(self.flux) * amp_factors[idx]
            self.lower_bounds[npoly + 3 * idx + 1] = 0
            self.upper_bounds[npoly + 3 * idx + 1] = np.max(self.flux)
            self.guess[npoly + 3 * idx + 2] = self.wave_c
            self.lower_bounds[npoly + 3 * idx + 2] = self.wave_c - 0.1
            self.upper_bounds[npoly + 3 * idx + 2] = self.wave_c + 0.1
            guess_width = self.wave[int(width_data[3])] - self.wave[int(width_data[2])] * broad_factors[idx]
            if guess_width == 0:
                self.guess[npoly + 3 * idx + 3] = 0.01
            else:
                self.guess[npoly + 3 * idx + 3] = guess_width * broad_factors[idx]
        
        #print(self.guess)
        if spaxel_fit:
            try: 
                popt, pcov = curve_fit(fitting_function, self.wave, self.flux, self.guess, bounds=(self.lower_bounds, self.upper_bounds), maxfev=5000)
                self.fitfunc = fitting_function
                
                self.popt = popt 
                self.pcov = pcov
                self.fitflag = "custom"
                self.residuals = (self.fitfunc(self.wave, *self.popt) - self.flux) / self.flux * 100
                
                #gamma = self.popt[4] * 2.355 / self.wave_c
                #print(gamma)
                #print(self.popt[2])
                #gauss = [popt[2], popt[3], gamma]
                #gauss = [popt[3], popt[4], popt[2]]
                #flux = gauss_prof(self.wave, gauss, ext=None)
                
                self.line_strength = 0
                for idx in range(ngauss):
                    gamma = np.abs(self.popt[npoly + 3 * idx + 3]) * 2.355 / self.popt[npoly + 3 * idx + 2]
                    #print(gamma, self.popt[npoly + 3 * idx + 1], self.popt[npoly + 3 * idx + 2], self.popt[npoly + 3 * idx + 3])
                    self.line_strength += (1 / np.sqrt(np.pi * np.log(2)) * (np.pi * const.c.to('micron/s') / 2) * (self.popt[npoly + 3 * idx + 1] * u.Jy * gamma / (self.popt[npoly + 3 * idx + 2] * u.micron))).to(u.watt / u.meter ** 2)
                        
                return popt, pcov
            except:
                self.popt = None 
                self.pcov = None
                self.line_strength = 0
                return None
        
        else:
            popt, pcov = curve_fit(fitting_function, self.wave, self.flux, self.guess, bounds=(self.lower_bounds, self.upper_bounds), maxfev=5000)
            self.fitfunc = fitting_function
            
            self.popt = popt 
            self.pcov = pcov
            self.fitflag = "custom"
            self.residuals = (self.fitfunc(self.wave, *self.popt) - self.flux) / self.flux * 100
            
            #gamma = self.popt[4] * 2.355 / self.wave_c
            #print(gamma)
            #print(self.popt[2])
            #gauss = [popt[2], popt[3], gamma]
            #gauss = [popt[3], popt[4], popt[2]]
            #flux = gauss_prof(self.wave, gauss, ext=None)
            
            self.line_strength = 0
            for idx in range(ngauss):
                gamma = np.abs(self.popt[npoly + 3 * idx + 3]) * 2.355 / self.popt[npoly + 3 * idx + 2]
                #print(gamma, self.popt[npoly + 3 * idx + 1], self.popt[npoly + 3 * idx + 2], self.popt[npoly + 3 * idx + 3])
                self.line_strength += (1 / np.sqrt(np.pi * np.log(2)) * (np.pi * const.c.to('micron/s') / 2) * (self.popt[npoly + 3 * idx + 1] * u.Jy * gamma / (self.popt[npoly + 3 * idx + 2] * u.micron))).to(u.watt / u.meter ** 2)
                    
            return popt, pcov
    
    
    def neon3_fit(self):
        """ 
        The main refitting routine. The continuum function is defined
        as an `n` order polynomial. The line is assumed to be 
        comprised of strictly gaussian components.
        
        Arguments
        ---------
        npoly : int
            The degree of the polynomial to be fit as the continuum 
            function. The default value is 1.
        ngauss : int 
            The number of gaussian components to be fit as the line 
            feature. The default value is 2.
        """
        
        ngauss = 3
        npoly = 1
        # Define the initial parameter map
        self.n_params = (npoly + 1) + ngauss * 3
        self.params = [None] * self.n_params
        self.param_map = {}
        self.param_map_r = {}
        self.guess = [0] * (npoly + ngauss * 3 + 1)
        self.lower_bounds = [-np.inf] * (npoly + ngauss * 3 + 1)
        self.upper_bounds = [np.inf] * (npoly + ngauss * 3 + 1)
        self.npoly = npoly 
        self.ngauss = ngauss
        
        # Construct the parameter map
        gauss_idx = 1
        ngauss_idx = 1
        for idx, _ in enumerate(self.params):
            if idx < (npoly + 1):
                self.param_map[f"poly_{npoly - idx}"] = idx
                self.param_map_r[str(idx)] = f"poly_{npoly - idx}"
            if idx >= (npoly + 1):
                self.param_map[f"gauss_{ngauss_idx}_{gauss_idx}"] = idx
                self.param_map_r[str(idx)] = f"gauss_{ngauss_idx}_{gauss_idx}"
                gauss_idx += 1
                if gauss_idx == 4:
                    gauss_idx = 1
                    ngauss_idx += 1
        
        # Define the continuum function
        if npoly == 1:
            self.continuum = OneDPolynomial
        if npoly == 2:
            self.continuum = TwoDPolynomial
        
        # Define the complete set of gaussians
        for _ in range(ngauss):
            self.gaussians.append(OneDGaussian)
        
        
        def fitting_function(wave, *pars):
            """ 
            Dynamically defined fitting function.
            """
            if self.npoly == 1:
                continuum = self.continuum(wave, pars[0], pars[1])
            if self.npoly == 2:
                continuum = self.continuum(wave, pars[0], pars[1], pars[2])
            gaussian_values = []
            for idx, gaussian in enumerate(self.gaussians):
                gaussian_values.append(gaussian(wave, pars[npoly + 3 * idx + 1], pars[npoly + 3 * idx + 2], pars[npoly + 3 * idx + 3]))
            fluxes = continuum
            for arrays in gaussian_values:
                fluxes += arrays
            return fluxes
        
        if npoly == 1:
            self.guess[0] = (self.flux[-1] - self.flux[0]) / (self.wave[-1] - self.wave[0])
            self.guess[1] = self.flux[0] - self.guess[0] * self.wave[0]
        if npoly == 2:
            self.guess[0] = 0
            self.guess[1] = 0
            self.guess[2] = 0
        
        amp_factors = [1, 1, 0.2, 0.1]
        broad_factors = [1, 1, 2, 3]
        width_data = peak_widths(self.flux, [self.peak_idx])

        wave_cs = [15.552, 15.555, 15.557]
        for idx in range(ngauss):
            self.guess[npoly + 3 * idx + 1] = np.max(self.flux) * amp_factors[idx]
            self.lower_bounds[npoly + 3 * idx + 1] = 0
            self.upper_bounds[npoly + 3 * idx + 1] = np.max(self.flux)
            self.guess[npoly + 3 * idx + 2] = self.wave_c
            self.lower_bounds[npoly + 3 * idx + 2] = wave_cs[idx] - 0.001
            self.upper_bounds[npoly + 3 * idx + 2] = wave_cs[idx] + 0.001
            guess_width = self.wave[int(width_data[3])] - self.wave[int(width_data[2])] * broad_factors[idx]
            if guess_width == 0:
                self.guess[npoly + 3 * idx + 3] = 0.01
            else:
                self.guess[npoly + 3 * idx + 3] = guess_width * broad_factors[idx]
        
        #print(self.guess)
        print(self.guess)
        if True:
            try: 
                popt, pcov = curve_fit(fitting_function, self.wave, self.flux, self.guess, bounds=(self.lower_bounds, self.upper_bounds), maxfev=5000)
                self.fitfunc = fitting_function
                
                self.popt = popt 
                self.pcov = pcov
                self.fitflag = "custom"
                self.residuals = (self.fitfunc(self.wave, *self.popt) - self.flux) / self.flux * 100
                
                #gamma = self.popt[4] * 2.355 / self.wave_c
                #print(gamma)
                #print(self.popt[2])
                #gauss = [popt[2], popt[3], gamma]
                #gauss = [popt[3], popt[4], popt[2]]
                #flux = gauss_prof(self.wave, gauss, ext=None)
                
                self.line_strength = 0
                for idx in range(ngauss):
                    gamma = np.abs(self.popt[npoly + 3 * idx + 3]) * 2.355 / self.popt[npoly + 3 * idx + 2]
                    #print(gamma, self.popt[npoly + 3 * idx + 1], self.popt[npoly + 3 * idx + 2], self.popt[npoly + 3 * idx + 3])
                    self.line_strength += (1 / np.sqrt(np.pi * np.log(2)) * (np.pi * const.c.to('micron/s') / 2) * (self.popt[npoly + 3 * idx + 1] * u.Jy * gamma / (self.popt[npoly + 3 * idx + 2] * u.micron))).to(u.watt / u.meter ** 2)
                        
                return popt, pcov
            except:
                self.popt = None 
                self.pcov = None
                self.line_strength = 0
                return None
        
        else:
            popt, pcov = curve_fit(fitting_function, self.wave, self.flux, self.guess, bounds=(self.lower_bounds, self.upper_bounds), maxfev=5000)
            self.fitfunc = fitting_function
            
            self.popt = popt 
            self.pcov = pcov
            self.fitflag = "custom"
            self.residuals = (self.fitfunc(self.wave, *self.popt) - self.flux) / self.flux * 100
            
            #gamma = self.popt[4] * 2.355 / self.wave_c
            #print(gamma)
            #print(self.popt[2])
            #gauss = [popt[2], popt[3], gamma]
            #gauss = [popt[3], popt[4], popt[2]]
            #flux = gauss_prof(self.wave, gauss, ext=None)
            
            self.line_strength = 0
            for idx in range(ngauss):
                gamma = np.abs(self.popt[npoly + 3 * idx + 3]) * 2.355 / self.popt[npoly + 3 * idx + 2]
                #print(gamma, self.popt[npoly + 3 * idx + 1], self.popt[npoly + 3 * idx + 2], self.popt[npoly + 3 * idx + 3])
                self.line_strength += (1 / np.sqrt(np.pi * np.log(2)) * (np.pi * const.c.to('micron/s') / 2) * (self.popt[npoly + 3 * idx + 1] * u.Jy * gamma / (self.popt[npoly + 3 * idx + 2] * u.micron))).to(u.watt / u.meter ** 2)
                    
            return popt, pcov
    
    
    def fit_test(self):
        """ 
        This functions runs a basic residual test to determine optimal fitting
        """
        if self.fitflag == None:
            raise AssertionError("Fit has not been performed!")
        
        abs_residuals = np.absolute(self.residuals)
        fig, ax = plt.subplots()
        ax.plot(self.wave, abs_residuals)
        plt.show()
    
    
    ### Visualization Functions
    
    
    def render_data(self, endaction="display"):
        """ 
        Renders the data presented to FitFuncs
        """
        
        #plt.style.use('dark_background')
        fig, ax = plt.subplots()
        fig.set_size_inches(8, 8)
        ax.scatter(self.wave, self.flux, s=10)
        ax.axvline(self.wave_c, ls="dotted", color="yellow", alpha=0.5)
        ax.set_xlim(self.wave_range[0], self.wave_range[1])
        ax.set_title(self.name, loc="right")
        ax.set_yscale("log")
        
        # Determine ending action for script
        if endaction == "display":
            plt.show()
        if endaction == "savefig":
            plt.savefig("./test.png")
        if endaction == "return":
            return fig, ax

    
    def render_fit(self, path):
        """ 
        Generates a two subplot diagram which displays the fitted 
        function, inputted data, and fourth standard deviation 
        residuals.
        """

        # Sanity check if fit was complete
        if self.fitflag == None:
            raise AssertionError("Fit has yet to be performed on inputted data!")
        
        # Define figure and axes
        fig, ((ax1), (ax2)) = plt.subplots(nrows=2,ncols=1, sharex=True, height_ratios=[3, 1])
        fig.subplots_adjust(hspace=0.0)
        
        # Plot data, total fit, and continuum
        ax1.scatter(self.wave, self.flux, color="gold", s=5, label="Spectrum")
        if self.npoly == 1:
            cont = OneDPolynomial(self.wave, self.popt[0], self.popt[1])
        if self.npoly == 2:
            cont = TwoDPolynomial(self.wave, self.popt[0], self.popt[1], self.popt[2])
        ax1.plot(self.wave, self.fitfunc(self.wave, *self.popt), color="white", ls="solid", linewidth=1, label='Fit')
        ax1.plot(self.wave, cont, color="white", ls="dashed", linewidth=1, label='Continuum')
        
        fit_width = np.max([self.popt[self.npoly + 3 * idx + 3] for idx in range(self.ngauss)])
        xlim_low = np.max([self.wave_c - fit_width * 7, self.wave[0]])
        xlim_high = np.min([self.wave_c + fit_width * 7, self.wave[-1]])
        ax1.set_xbound(xlim_low, xlim_high)
        # Plot gaussians
        for idx, _ in enumerate(self.gaussians):
            ax1.plot(self.wave, OneDGaussian(self.wave, self.popt[self.npoly + 3 * idx + 1], self.popt[self.npoly + 3 * idx + 2], self.popt[self.npoly + 3 * idx + 3]) + cont, ls=self.linestyles[idx], color=self.plotcolors[idx], linewidth=1, label=f'Gaussian {idx + 1}')
        
        # Plot residuals
        ax2.plot(self.wave, self.residuals, color="white")
        
        # Set axes labels
        ax2.set_xlabel(r"$\lambda$ (microns)")
        ax1.set_ylabel(r"$f_{\nu} (Jy)$")
        ax2.set_ylabel("Residuals (%)")
        
        # Define plot limits for visibility
        std = np.nanstd(self.residuals)
        
        
        
        ax1.set_ylim(np.min(self.flux) * 0.8, np.max(self.flux) * 1.2)
        ax1.set_yscale("log")
        
        ax2.set_ylim(-4*std, 4*std)
        ax1.legend()
        
        #plt.show()
        plt.savefig(path, dpi=1000, bbox_inches="tight")
        plt.close()