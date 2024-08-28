import os 
import glob
import site
from textwrap import dedent
import numbers

import numpy as np
import astropy.io.fits as fits
import astropy.units as u
import astropy.constants as const
from astropy import coordinates
from astroquery.ipac.ned import Ned
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt
import matplotlib.animation as animation


import CRETA.creta as creta
import creta_mod 

import CAFE.cafe as cafe

from localfit import LocalFit
from plotparams import PlotParams

import asdf
from asdf import AsdfFile
from CAFE.component_model import pah_drude, gauss_prof, drude_prof, drude_int_fluxes


class CubeSpec:
    """ 
    A container for performing spectral extraction and fitting 
    operations with CRETA and CAFE.
    """
    
    def __init__(self, creta_dir, param_path, param_file, creta_input_path, 
                 redshift=None, fit_dirname=None, mode="SB_and_AGN"):
        
        # Set and verify essential path variables
        self._creta_dir = os.path.abspath(creta_dir)
        self.creta_input_path = os.path.join(self._creta_dir, 
                                             creta_input_path)
        self._param_path = os.path.join(self._creta_dir, param_path)
        self.fit_dirname = fit_dirname
        
        
        self._verify_paths()
        
        # Set and verify essential files
        self.param_file = param_file
        self._verify_files()
        
        # Set essential class fields
        self.c = None
        self.redshift = redshift
        self.mode = mode
        self.c_mod = None
        
        
        # Query for a redshift using NED if none was provided
        if self.redshift == None:
            self.find_redshift()
    
    
    ### Housekeeping methods
    
    
    def generate_params(self):
        pass
    
    
    def _verify_paths(self):
        """ 
        Checks if all inputted paths are present prior to extraction
        or fitting and either generates such directories or raises an
        error if a directory containing apriori data is missing.
        """
        
        if not os.path.exists(self._creta_dir):
            error_string = f"""
            Specified CRETA directory:
            {self._creta_dir}
            does not exist or is not correctly inputted!
            """
            raise NotADirectoryError(dedent(error_string))
        if not os.path.exists(self.creta_input_path):
            error_string = f"""
            Specified CRETA input directory:
            {self.creta_input_path}
            does not exist or is not correctly inputted!
            """
            raise NotADirectoryError(dedent(error_string))
        
        if not os.path.exists(self._param_path):
            error_string = f"""
            Specified parameter file directory:
            {self._param_path}
            does not exist or is not correctly inputted!
            """
            raise NotADirectoryError(dedent(error_string))
    
    
    def _verify_files(self):
        """ 
        Verifies that inputted data files are present for running this 
        pipeline. This method assumes that only the Level3-reduced 
        datacubes are present in the input data folder and that the 
        parameter file files a certain naming structure.
        """
        
        # List all data files
        self.data_files = glob.glob(os.path.join(self.creta_input_path,
                                                 "*.fits"))
        test_header = fits.getheader(self.data_files[0], ext=0)
        test_target = test_header["TARGPROP"]
        
        # Check parameter files
        with open(os.path.join(self._param_path, self.param_file)) as file:
            param_lines = file.readlines()
        self.param_dict = {}
        if "single" in self.param_file:
            self.param_dict["type"] = "single"
        if "grid" in self.param_file:
            self.param_dict["type"] = "grid"
        self.extraction_id = self.param_file.split("_")[1]
        for idx, line in enumerate(param_lines):
            if idx == 0:
                key = line.split("=")[0]
                values = line.split("=")[1].split(",")
                cubes = []
                for value in values:
                    cubes.append(value.strip())
                self.param_dict[key] = cubes
            if idx != 0:
                key = line.split("=")[0].strip()
                value = line.split("=")[1].split("#")[0].strip()
                if key == "user_r_ap":
                    value = float(value)
                    self.asec = value
                if key == "lambda_ap":
                    value = float(value)
                if key == "lambda_cent":
                    value = float(value)
                if key == "r_ann_in":
                    value = float(value)
                if key == "ann_width":
                    value = float(value)
                self.param_dict[key] = value
        #print(self.param_dict)
        
        # Cross check all data files
        if len(self.data_files) == 12:
            for filename in self.data_files:
                current_header = fits.getheader(filename, ext=0)
                current_target = current_header["TARGPROP"]
                if current_target != test_target:
                    raise ValueError("Datacubes are of different targets!")
            self.ra = test_header["TARG_RA"]
            self.dec = test_header["TARG_DEC"]
            self.target = test_target
            self.header = test_header
            
            # Verify existence of parameter file
            param_file_path = os.path.join(self._param_path, self.param_file)
            if not os.path.exists(param_file_path):
                error_string = f"""
                Specified parameter file:
                {param_file_path}
                does not exist or is not correctly inputted!
                """
                raise FileNotFoundError(dedent(error_string))
            
            # Verify existance of CRETA and CAFE output path
            creta_output_extension = f"creta_output/extractions/{self.target}/{self.extraction_id}"
            self.creta_output_path = os.path.join(self._creta_dir, 
                                              creta_output_extension)
            if not os.path.exists(self.creta_output_path):
                os.makedirs(self.creta_output_path)
                method_string = f""" 
                Specified CRETA output directory:
                {self.creta_output_path}
                did not exist and has been generated.
                """
                print(method_string)
            cafe_output_extension = f"cafe_output/{self.target}/{self.extraction_id}/{self.fit_dirname}"
            self.cafe_output_path = os.path.join(self._creta_dir,
                                                 cafe_output_extension)
            if not os.path.exists(self.cafe_output_path):
                os.makedirs(self.cafe_output_path)
                method_string = f""" 
                Specified CAFE output directory:
                {self.cafe_output_path}
                did not exist and has been generated.
                """
                print(method_string)
        else:
            if len(self.data_files) < 12:
                raise FileNotFoundError("Datacubes missing!")
            if len(self.data_files) > 12:
                raise FileExistsError("Extraneous files found!")
    
    
    def _initialize_CRETA(self):
        """ 
        Initializes an instance of CRETA for this particular MIRI 
        observation.
        """
        self.c = creta.creta(self._creta_dir)
        return True
    
    
    # EXPERIMENTAL METHOD
    
    
    def _initialize_CRETA_mod(self):
        """ 
        Initializes a modified instance of CRETA for this particular MIRI 
        observation.
        """
        self.c_mod = creta_mod.creta(self._creta_dir)
        return True
    
    
    def rewrite_spec_csv(self):
        """ 
        This method rewrites the CRETA outputted csv into one that is
        accepted by Thomas Lai's web app for viewing spectra.
        """
        import csv 
        
        method_string = f"""
        Reformating spectrum csv file to be Thomas Lai compliant
        """
        print(dedent(method_string))
        
        new_asec = str(self.asec)[0] + "." + str(self.asec)[1]
        original_csv = os.path.join(self.creta_output_path, 
                                    f"{self.target}_SingleExt_r{new_asec}as.csv")
        spec_dict = {"w": [], "f": [], "f_unc": []}
        readlines = False
        stop_string = "Wave,Band_name,Flux_ap,Err_ap,R_ap,Flux_ap_st,"
        
        with open(original_csv, 'r') as csvfile:
            for line in csvfile.readlines():
                # Ignore up to this line
                if stop_string in line:
                    readlines = True 
                    continue
                # Record values
                if readlines:
                    vals = line.split(sep=",")
                    spec_dict["w"].append(vals[0])
                    spec_dict["f"].append(vals[5])
                    spec_dict["f_unc"].append(vals[6])

        new_csv = os.path.join(self.creta_output_path, 
                               f"{self.target}_sum_sf_spec.csv")

        with open(new_csv, 'w') as csvfile:
            fieldnames = ["w", "f", "f_unc"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for idx, _ in enumerate(spec_dict["w"]):
                writer.writerow({"w": spec_dict["w"][idx], 
                                 "f" : spec_dict["f"][idx], 
                                 "f_unc" : spec_dict["f_unc"][idx]})

        method_string = f"""
        New csv file found at: {new_csv}
        """
        print(dedent(method_string))
        return None
    
    
    def find_redshift(self, query_rad=5):
        """ 
        Redshift autoquery routine using the NED astroquery interface.
        
        Arguments
        ---------
        query_rad : float
            The size of the query circle in arcseconds.
        """
        
        query_coord = coordinates.SkyCoord(ra=self.header["TARG_RA"],
                                           dec=self.header["TARG_DEC"],
                                           unit=(u.deg, u.deg))
        result_table = Ned.query_region(query_coord, 
                                        radius=query_rad*u.arcsec)
        result_table = result_table[~np.isnan(result_table['Redshift'])]   
        n_results = len(result_table)
        
        if n_results == 1:
            self.redshift = result_table["Redshift"][0]
            method_string = f"""
            Autoquerying for redshift using header!
            RA={self.ra}, DEC={self.dec}
            Found {n_results} result(s): z={self.redshift}
            """
            print(dedent(method_string))
            return self.redshift
        if n_results == 0:
            error_string = f"""
            Redshift not found. Try increasing the query radius 
            parameter to get a result.
            """
            raise ValueError(dedent(error_string))
        if n_results > 1:
            error_string = f"""
            Multiple sources found within query radius. Try
            decreasing with the query radius parameter to
            reduce returned results.
            """
            print(result_table)
            raise ValueError(dedent(error_string))
    
    
    ### Value recovery routines
    
    
    def open_asdf(self):
        """ 
        This routine opens the asdf file corresponding to the last 
        CAFE fitting session.
        """
        asdf_fn = os.path.join(self.cafe_output_path, f"{self.target}_SingleExt_r{str(self.asec).replace('.', '')}as/{self.target}_SingleExt_r{str(self.asec).replace('.', '')}as_cafefit.asdf")
        print("ASDF", asdf_fn)
        af = asdf.open(asdf_fn)
        return af
    
    
    def recall_data(self):
        """ 
        Recalls the saved spectral profile parameters for a previous
        CAFE fitting session.
        """
        
        af = self.open_asdf()
        wave = np.asarray(af.tree['cafefit']['obsspec']['wave'])
        flux = np.asarray(af['cafefit']['obsspec']['flux'])
        flux_unc = np.asarray(af['cafefit']['obsspec']['flux_unc'])
        
        spec_dict = {'wave':wave, 'flux':flux, 'flux_unc':flux_unc}
        return spec_dict
    
    
    def recall_comps(self):
        """ 
        Recalls the saved component parameters for a previous CAFE
        fitting session.
        """
        
        af = self.open_asdf()
        comps = af['cafefit']['CompFluxes']
        for key in comps.keys():
            comps[key] = np.asarray(comps[key])
        return comps
    
    
    def recall_gauss(self):
        """ 
        Recalls the saved gaussian profile parameters for a previous
        CAFE fitting session.
        """
        
        af = self.open_asdf()
        g = af['cafefit']['gauss']
        gauss = [np.asarray(g['wave']), np.asarray(g['gamma']), np.asarray(g['peak']), np.asarray(g['name'])]
        return gauss
    
    
    def recall_drude(self):
        """ 
        Recalls the saved drude profile parameters for a previous
        CAFE fitting session.
        """
        
        af = self.open_asdf()
        d = af['cafefit']['drude']
        drude = [np.asarray(d['wave']), np.asarray(d['gamma']), np.asarray(d['peak']), np.asarray(d['name'])]
        return drude
    
    
    def recall_extPAH(self):
        """ 
        Recalls the saved PAH extinction parameters for a previous CAFE
        fitting session.
        """
        
        af = self.open_asdf()
        extPAH = np.asarray(af['cafefit']['extComps']['extPAH'])
        return extPAH
    
    
    ### Analysis routines
    
    
    def perform_single_extraction(self):
        """ 
        This method runs the CRETA single extraction tool on this MIRI
        observation. Assumes an existing single parameter file exists.
        """
        if self.c == None:
            self._initialize_CRETA()
        if self.param_dict["type"] != "single":
            raise ValueError("Parameter file is not configured for single extraction.")
        self.c.singleExtraction(parameter_file=True, 
                                parfile_path=self._param_path, 
                                parfile_name="/" + self.param_file, 
                                data_path=self.creta_input_path, 
                                output_path=self.creta_output_path + "/",
                                output_filebase_name=f'{self.target}',
                                ignore_DQ=True)
    
    
    # EXPERIMENTAL METHOD
    
    def perform_single_extraction_custom(self):
        """ 
        This method runs a modified CRETA single extraction tool on this MIRI
        observation. Assumes an existing single parameter file exists.
        """
        if self.c_mod == None:
            self._initialize_CRETA_mod()
        if self.param_dict["type"] != "single":
            raise ValueError("Parameter file is not configured for single extraction.")
        realData = self.c_mod.singleExtraction(parameter_file=True, 
                                parfile_path=self._param_path, 
                                parfile_name="/" + self.param_file, 
                                data_path=self.creta_input_path, 
                                output_path=self.creta_output_path + "/",
                                output_filebase_name=f'{self.target}',
                                ignore_DQ=True)
        return realData
    
    
    def perform_grid_extraction(self):
        """ 
        This method runs the CRETA grid extraction tool on this MIRI
        observation. Assumes an existing grid parameter file exists.
        """
        if self.c_mod == None:
            self._initialize_CRETA_mod()
        if self.param_dict["type"] != "grid":
            raise ValueError("Parameter file is not configured for single extraction.")
        self.c_mod.gridExtraction(data_path=self.creta_input_path, 
                              parfile_path=self._param_path, 
                              parfile_name="/" + self.param_file, 
                              output_path = self.creta_output_path + "/",
                              output_filebase_name=f'{self.target}',
                              ignore_DQ=True)
    
    
    
    
    
    def perform_fit(self):
        """ 
        This method runs the CAFE fitting tool on this MIRI observation
        given that CRETA outputs are already available for this object.
        """
        
        # Define appropriate directories and pathing structures
        cafe_site_dir = site.getsitepackages()[0] + "/CAFE/"
        print(self.target)
        cafe_output_path = os.path.join(self._creta_dir, f"cafe_output/{self.target}/{self.extraction_id}/{self.fit_dirname}")
        source_fd = self.creta_output_path
        new_asec = str(self.asec)[0] + "." + str(self.asec)[1]
        #print(new_asec)
        #print(self.asec)
        source_fn = f'{self.target}_SingleExt_r{self.asec}as.fits'
        #inppar_fn = cafe_site_dir + "inp_parfiles/inpars_jwst_miri_AGN.ini"
        #optpar_fn = cafe_site_dir + "opt_parfiles/default_opt.cafe"
        
        inppar_fn = cafe_site_dir + f"inp_parfiles/inpars_jwst_miri_{self.mode}.ini"
        optpar_fn = cafe_site_dir + "opt_parfiles/default_opt.cafe"
        inppar_fn = os.path.join(self._creta_dir, f"inp_parfiles/inpars_jwst_miri_{self.mode}.ini")
        optpar_fn = os.path.join(self._creta_dir, "opt_parfiles/default_opt.cafe")
        print(inppar_fn)
        print(optpar_fn)
        
        print(cafe_output_path)
        # Initialize the spectroscopic fitting class
        s = cafe.specmod(cafe_site_dir)
        
        # Read in the spectrum
        s.read_spec(source_fn, file_dir=source_fd + '/', z=self.redshift)

        # Preview the spectrum and overlay the initial params
        #s.plot_spec_ini(inppar_fn, optpar_fn)
        
        # Thoroughly fit the spectrum
        s.fit_spec(inppar_fn, optpar_fn, output_path=cafe_output_path + '/')
        
        # Plot the spectrum
        #s.plot_spec_fit(inppar_fn, optpar_fn)
        #self.cafeplot()
    
    
    def return_lines(self):
        return os.path.join(self.cafe_output_path, f"{self.target}_SingleExt_r{str(self.asec).replace('.', '')}as/{self.target}_SingleExt_r{str(self.asec).replace('.', '')}as_linetable.ecsv")
    
    def return_pah(self):
        return os.path.join(self.cafe_output_path, f"{self.target}_SingleExt_r{str(self.asec).replace('.', '')}as/{self.target}_SingleExt_r{str(self.asec).replace('.', '')}as_pahtable.ecsv")
    
    def run_spectool(self):
        """ 
        This script runs Thomas Lai's spectool and can be used to make
        bird's eye observations of a given spectra and any stand-out 
        line features.
        """
        import subprocess 
        
        os.chdir("./../JWST-SpecTool/src/")
        
        subprocess.run(["python", "-m", "webbrowser", "-t", "http://127.0.0.1:8050"])
        subprocess.run(["python", "app.py", f"{self.redshift}", f"{self.target}"])
    
    
    def recall_fit(self):
        """ 
        Recalls the saved fit parameters for a previous CAFE fitting 
        session.
        """
        
        asdf_fn = os.path.join(self.cafe_output_path, f"{self.target}_SingleExt_r{str(self.asec).replace('.', '')}as/{self.target}_SingleExt_r{str(self.asec).replace('.', '')}as_cafefit.asdf")
        print(asdf_fn)
        af = asdf.open(asdf_fn)
        
        wave = np.asarray(af.tree['cafefit']['obsspec']['wave'])
        flux = np.asarray(af['cafefit']['obsspec']['flux'])
        flux_unc = np.asarray(af['cafefit']['obsspec']['flux_unc'])
        
        comps = af['cafefit']['CompFluxes']
        extPAH = af['cafefit']['extComps']['extPAH']
        g = af['cafefit']['gauss']
        d = af['cafefit']['drude']
        
        gauss = [g['wave'], g['gamma'], g['peak']]
        drude = [d['wave'], d['gamma'], d['peak']]
        
        spec_dict = {'wave':wave, 'flux':flux, 'flux_unc':flux_unc}
        
        # Assuming there is no phot input.
        # TODO: include phot_dict as input.
        self.cafeplot(spec_dict, None, comps, gauss, drude, pahext=extPAH)
        pass
    
    
    def recall_line(self, low_lamb=None, high_lamb=None):
        """ 
        Recalls the saved fit parameters for a previous CAFE fitting 
        session.
        """
        
        asdf_fn = os.path.join(self.cafe_output_path, f"{self.target}_SingleExt_r{str(self.asec).replace('.', '')}as/{self.target}_SingleExt_r{str(self.asec).replace('.', '')}as_cafefit.asdf")
        print(asdf_fn)
        af = asdf.open(asdf_fn)
        
        wave = np.asarray(af.tree['cafefit']['obsspec']['wave'])
        flux = np.asarray(af['cafefit']['obsspec']['flux'])
        flux_unc = np.asarray(af['cafefit']['obsspec']['flux_unc'])
        
        comps = af['cafefit']['CompFluxes']
        extPAH = af['cafefit']['extComps']['extPAH']
        g = af['cafefit']['gauss']
        d = af['cafefit']['drude']
        
        gauss = [g['wave'], g['gamma'], g['peak']]
        drude = [d['wave'], d['gamma'], d['peak']]
        
        spec_dict = {'wave':wave, 'flux':flux, 'flux_unc':flux_unc}
        line_names = af['cafefit']['gauss']['name']
        line_waves = np.asarray(af['cafefit']['gauss']['wave'])
        print(af['cafefit']['gauss'].keys())
        
        # Assuming there is no phot input.
        # TODO: include phot_dict as input.
        self.line_diagnostic(spec_dict, None, comps, gauss, drude, line_names, line_waves, low_lamb, high_lamb, pahext=extPAH)
        #pass
    
    
    def local_fit_flux(self, wave_range, wave_c, name, npoly, ngauss):
        """ 
        This routine generates a local flux estimate for a specified line.
        """
        data = self.recall_data()
        figob = LocalFit(data["wave"], data["flux"], wave_range, wave_c, name)
        figob.main_fit(npoly, ngauss)
        path = os.path.join(self.cafe_output_path, f"refit_{name}_{npoly}_{ngauss}.pdf")
        figob.render_fit(path)
        return figob.line_strength
    
    
    def cafeplot(self, spec, phot, comps, gauss, drude, vgrad={'VGRAD':0.}, plot_drude=True, pahext=None, save_name=False, params=None):
        ''' Plot the SED and the CAFE fit over the spectrum wavelength range

        Arguments:
        wave -- rest wavelength of observed spectrum
        flux -- observed flux values
        func -- uncertainties in measured fluxes
        comps -- dict of component fluxes

        Keyword Arugments:
        weights -- CAFE weights to use in estimating final chi^2 (default None)
        drude -- The collection of ouput parameters of Drude profiles
        plot_drude -- if true, plots individual drude profiles, otherwise plots total
        PAH contribution. (default false)
        pahext -- if not None, applies extinction curve to PAHs

        Returns: Figure
        '''
        
        plt.style.use('dark_background')
        fCir = comps['fCIR']
        fCld = comps['fCLD']
        fCoo = comps['fCOO']
        fWrm = comps['fWRM']
        fHot = comps['fHOT']
        fStb = comps['fSTB']
        fStr = comps['fSTR']
        fDsk = comps['fDSK']
        fLin = comps['fLIN']
        fPAH = comps['fPAH']
        fMod = fCir + fCld + fCoo + fWrm + fHot + fStb + fStr + fDsk + fLin + fPAH
        fCont = fCir + fCld + fCoo + fWrm + fHot + fStb + fStr + fDsk
        wavemod = comps['wave']
        
        fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios':[3,1]}, figsize=(8,8), sharex=True)
        ax1.scatter(spec['wave'], spec['flux'], color="white", s=2, edgecolor='white', facecolor='none', label='Spec Data', alpha=0.6, zorder=0)
        #ax1.errorbar(spec['wave'], spec['flux'], yerr=spec['flux_unc'], fmt='none', color='white', alpha=0.1)
        if phot is not None:
            ax1.scatter(phot['wave'], phot['flux'], marker='x', s=18, edgecolor='none', facecolor='black', label='Phot Data', alpha=0.9)
            ax1.errorbar(phot['wave'], phot['flux'], xerr=phot['width']/2, yerr=phot['flux_unc'], fmt='none', color='white', alpha=0.1)
            wave = np.concatenate((spec['wave'], phot['wave']))
            flux = np.concatenate((spec['flux'], phot['flux']))
            sortinds = np.argsort(wave)
            wave = wave[sortinds] ; flux = flux[sortinds]
        else:
            wave = spec['wave']
            flux = spec['flux']
                                
        ax1.plot(wavemod, fCont, color='white', label='Continuum Fit', linestyle="dashed", zorder=4, lw=2, alpha=0.8)
        ax1.plot(wavemod, fCont+fLin+fPAH, color='cyan', label='Total Fit', linewidth=1.5, zorder=5, alpha=1) # green"""

        CLD_TMP = '' if params == None else r' ('+"{:.0f}".format(params['CLD_TMP'].value)+'$\,$K'+')'
        COO_TMP = '' if params == None else r' ('+"{:.0f}".format(params['COO_TMP'].value)+'$\,$K'+')'
        WRM_TMP = '' if params == None else r' ('+"{:.0f}".format(params['WRM_TMP'].value)+'$\,$K'+')'
        HOT_TMP = '' if params == None else r' ('+"{:.0f}".format(params['HOT_TMP'].value)+'$\,$K'+')'
            
        alpha = 1
        lw = 2
        if np.any(fCir > 0):
            ax1.plot(wavemod, fCir, label='Cirrus', c='tab:cyan', alpha=alpha, linewidth=lw)
        if np.sum(fCld > 0):
            ax1.plot(wavemod, fCld, label='Cold'+CLD_TMP, c='tab:blue', alpha=alpha, linewidth=lw, linestyle="dashed")
        if np.any(fCoo > 0):
            ax1.plot(wavemod, fCoo, label='Cool'+COO_TMP, c='green', alpha=alpha, linewidth=lw, linestyle="dotted") # teal
        if np.any(fWrm > 0):
            ax1.plot(wavemod, fWrm, label='Warm'+WRM_TMP, c='tab:orange', alpha=alpha, linewidth=lw, linestyle="dashdot")
        if np.any(fHot > 0):
            ax1.plot(wavemod, fHot, label='Hot'+HOT_TMP, c='brown', alpha=alpha, linewidth=lw, linestyle="dashed") # gold
        if np.any(fStb > 0): 
            ax1.plot(wavemod, fStb, label='Starburst', c='#FFEC00', alpha=alpha, linewidth=lw, linestyle="dotted")
        if np.any(fStr > 0):
            ax1.plot(wavemod, fStr, label='Stellar', c='#FF4500', alpha=alpha, linewidth=lw) # orangered
        #if np.any(fDsk > 0):
        #    ax1.plot(wavemod, fDsk, label='AGN', c='tab:red', alpha=alpha, linewidth=lw)
        if np.any(fLin > 0):
            ax1.plot(wavemod, fCont+fLin, label='Lines', c='lime', alpha=alpha, linewidth=lw, zorder=0) # blue ##1e6091

        """# Plot lines
        for i in range(len(gauss[0])):
            if pahext is None:
                pahext = np.ones(wavemod.shape)
            #print(gauss[0])
            lflux = gauss_prof(wavemod, [[gauss[0][i]], [gauss[1][i]], [gauss[2][i]]], ext=pahext)
            
            ax1.plot(wavemod, lflux+fCont, color='lime', label='_nolegend_', alpha=alpha, linewidth=1)
            if i == 0:
                ax1.plot(wavemod, lflux+fCont, color='lime', label='Lines', alpha=alpha, linewidth=1)
            #else:
                ax1.plot(wavemod, lflux+fCont, color='lime', label='_nolegend_', alpha=alpha, linewidth=1)
        """
        # Plot PAH features
        if plot_drude is True:
            for i in range(len(drude[0])):
                if pahext is None: 
                    pahext = np.ones(wavemod.shape)
                dflux = drude_prof(wavemod, [[drude[0][i]], [drude[1][i]], [drude[2][i]]], ext=pahext)

                if i == 0:
                    ax1.plot(wavemod, dflux+fCont, color='fuchsia', label='PAHs', alpha=alpha, linewidth=1)
                else:
                    ax1.plot(wavemod, dflux+fCont, color='fuchsia', label='_nolegend_', alpha=alpha, linewidth=1)
        elif np.any(fPAH > 0):
            ax1.plot(wavemod, fCont+fPAH, label='PAHs', color='fuchsia', alpha=alpha)

        #ax11 = ax1.twinx()
        #ax11.plot(wavemod, pahext, linestyle='dashed', color='pink', alpha=1, linewidth=1, label="Attenuation")
        #ax11.set_ylim(0, 1.1)
        #ax11.set_ylabel(r'Attenuation fraction $_{\rm{Warm\,dust, PAHs, Lines}}$', fontsize=14)
        #ax11.tick_params(axis='y', labelsize=10)
        #ax11.tick_params(direction='in', which='both', length=4, width=0.8, right=True)

        min_flux = np.nanmin(spec['flux'][np.r_[0:5,-5:len(spec['flux'])]])
        max_flux = np.nanmax(spec['flux'][np.r_[0:5,-5:len(spec['flux'])]])

        
        ax1.legend(loc='lower right')
        ax1.tick_params(direction='in', which='both', length=6, width=1, top=True)
        ax1.tick_params(axis='x', labelsize=0)
        ax1.tick_params(axis='y', labelsize=18)
        ax1.set_ylim(bottom=0.1*np.nanmin(min_flux), top=2.*np.nanmax(max_flux))
        ax1.set_xlim(left=4.7, right=15)
        ax1.set_xlim(np.nanmin(wave)/1.2, 1.2*np.nanmax(wave))
        ax1.set_ylabel(r'$f_\nu$ (Jy)', fontsize=20)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        #ax1.axvline(9.7, linestyle='--', alpha=0.2)

        xlabs = [1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50, 100, 200, 500]
        #ax1.set_xticks(xlabs[(np.where(xlabs > np.nanmin(wave))[0][0]):(np.where(xlabs < np.nanmax(wave))[0][-1]+1)])
        #ax1.xaxis.set_major_formatter(ScalarFormatter())

        interpMod = np.interp(wave, comps['wave'], fMod)
        res = (flux-interpMod) / flux * 100 # in percentage
        std = np.nanstd(res)
        ax2.plot(wave, res, color='white', linewidth=1)
        #ax2.plot(wave, (spec['flux']-interpMod)/func, color='k')
        ax2.axhline(0., color='black', linestyle='--')
        ax2.tick_params(direction='in', which='both', length=6, width=1,  right=True, top=True)
        ax2.tick_params(axis='x', labelsize=16)
        ax2.tick_params(axis='y', labelsize=16)
        ax2.set_ylim(-4*std, 4*std)
        #ax2.set_ylim(bottom=-4, top=4)
        ax2.set_xlabel(r'$\lambda_{\rm{rest}}$ $(\mu \rm{m})$', fontsize=20)
        #ax2.set_ylabel(r'$f^{data}_\nu - f^{tot}_\nu$ $(\sigma)$', fontsize=14)
        ax2.set_ylabel('Residuals (%)', fontsize=20)
        #ax1.set_zorder(100)
        ax1.set_xlim(4.7, 15)
        ax1.set_ylim(0.005, 0.3)

        fig.set_size_inches(12, 8)
        ax1.set_title(f'Southern Nucleus (0.3 asec)', loc="right", fontsize=16)
        plt.subplots_adjust(hspace=0)
        
        
        
        if save_name is True:
            plt.show()
            return (fig, ax1, ax2)
        else:
            #plot_path = os.path.join(self.cafe_output_path, f"{self.target}_SingleExt_r{self.asec}as/{self.target}_SingleExt_r{self.asec}as_fit5.pdf")
            plot_path = "pres_plot.pdf"
            fig.savefig(plot_path, dpi=1200, format='pdf', bbox_inches='tight')
            plt.close()
    
    
    def line_diagnostic(self, spec, phot, comps, gauss, drude, line_names, line_waves, lamb_low=None, lamb_high=None, vgrad={'VGRAD':0.}, plot_drude=True, pahext=None, save_name=False, params=None ):
        """ 
        Plots the SED and all fitted emission lines for a given target
        within some wavelength range.
        
        Arguments
        ---------
        lamb_low : float
            The lower limit of wavelength to be plotted.
        lamb_high : float
            The upper limit of wavelength to be plotted.
        """
        
        # Set range parameter and flag for non numeric values
        if lamb_low == None:
            lamb_low = min(spec['wave'] - 0.1)
        if lamb_high == None:
            lamb_high = max(spec['wave'] + 0.1)
        if not isinstance(lamb_low, numbers.Number):
            raise ValueError("Provided lower wavelength limit is not a number!")
        if not isinstance(lamb_high, numbers.Number):
            raise ValueError("Provided upper wavelength limit is not a number!")
        
        # Aggregate all continuum fit values 
        plt.style.use('dark_background')
        fCir = comps['fCIR']
        fCld = comps['fCLD']
        fCoo = comps['fCOO']
        fWrm = comps['fWRM']
        fHot = comps['fHOT']
        fStb = comps['fSTB']
        fStr = comps['fSTR']
        fDsk = comps['fDSK']
        fLin = comps['fLIN']
        fPAH = comps['fPAH']
        fMod = fCir + fCld + fCoo + fWrm + fHot + fStb + fStr + fDsk + fLin + fPAH
        fCont = fCir + fCld + fCoo + fWrm + fHot + fStb + fStr + fDsk
        wavemod = comps['wave']
        
        # Instantiate plot
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, 
                                       gridspec_kw={'height_ratios':[3,1]}, 
                                       figsize=(8,8), 
                                       sharex=True)
        ax1.scatter(spec['wave'], spec['flux'], color="white", s=2, 
                    edgecolor='white', facecolor='none', label='Spec Data', 
                    alpha=0.9, zorder=0)
        ax1.errorbar(spec['wave'], spec['flux'], yerr=spec['flux_unc'], 
                     capsize=2, capthick=1, color='white', 
                     alpha=1)
        if phot is not None:
            ax1.scatter(phot['wave'], phot['flux'], marker='x', s=18, edgecolor='none', facecolor='white', label='Phot Data', alpha=0.9)
            ax1.errorbar(phot['wave'], phot['flux'], xerr=phot['width']/2, yerr=phot['flux_unc'], fmt='none', color='white', alpha=0.1)
            wave = np.concatenate((spec['wave'], phot['wave']))
            flux = np.concatenate((spec['flux'], phot['flux']))
            sortinds = np.argsort(wave)
            wave = wave[sortinds] ; flux = flux[sortinds]
        else:
            wave = spec['wave']
            flux = spec['flux']

        fTot = fCont+fLin+fPAH
        ax1.plot(wavemod, fCont, color='white', label='Continuum Fit', linestyle="dashed", zorder=2, alpha=0.8)
        ax1.plot(wavemod, fTot, color='cyan', label='Total Fit', linewidth=1.5, zorder=5, alpha=0.85) # green
        

        # Plot lines
        for i in range(len(gauss[0])):
            if pahext is None:
                pahext = np.ones(wavemod.shape)

            lflux = gauss_prof(wavemod, [[gauss[0][i]], [gauss[1][i]], [gauss[2][i]]], ext=pahext)
            
            if i == 0:
                ax1.plot(wavemod, lflux+fCont, color='red', label='Lines', alpha=1, linewidth=2, zorder=4)
            else:
                ax1.plot(wavemod, lflux+fCont, color='red', label='_nolegend_', alpha=1, linewidth=2, zorder=4)
        
        skip_flag = False
        
        # Label lines
        for idx, wavenumber in enumerate(line_waves):
            if skip_flag:
                skip_flag = False
                continue
            if idx != len(line_waves) - 1:
                #print(line_names[idx], line_names[idx + 1])
                if line_names[idx] == line_names[idx + 1]:
                    #print("skipped line", line_names[idx])
                    skip_flag = True
            #if idx == len(line_waves):
                #print(line_names[idx])
            difference_array = np.absolute(wavemod - wavenumber)
            index = difference_array.argmin()
            ax1.axvline(wavenumber, color="yellow", zorder=0, alpha=0.5)
            name = line_names[idx].split('_')[0]
            if "H2" in name:
                formatted_name = r"H$_2$ 0-0 $S$" + f"({name[-1:]})"
                ax1.text(wavemod[index], fTot[index] * 1.1, formatted_name, bbox=dict(boxstyle="round",
                   ec=(1., 1, 1),
                   fc=(0., 0, 1),
                   ))
            else:
                formatted_name = name
                ax1.text(wavemod[index], fTot[index] * 1.1, formatted_name, bbox=dict(boxstyle="round",
                    ec=(1., 1, 1),
                    fc=(0., 0, 1),
                    ))
        
        # Plot PAH features
        if plot_drude is True:
            for i in range(len(drude[0])):
                if pahext is None: 
                    pahext = np.ones(wavemod.shape)
                dflux = drude_prof(wavemod, [[drude[0][i]], [drude[1][i]], [drude[2][i]]], ext=pahext)

                if i == 0:
                    ax1.plot(wavemod, dflux+fCont, color='pink', label='PAHs', alpha=1, linewidth=1, zorder=3)
                else:
                    ax1.plot(wavemod, dflux+fCont, color='pink', label='_nolegend_', alpha=1, linewidth=1, zorder=3)
        elif np.any(fPAH > 0):
            ax1.plot(wavemod, fCont+fPAH, label='PAHs', color='pink', alpha=1, zorder=3)

        min_flux = np.nanmin(spec['flux'][np.r_[0:5,-5:len(spec['flux'])]])
        max_flux = np.nanmax(spec['flux'][np.r_[0:5,-5:len(spec['flux'])]])

        ax1.legend(loc='lower right')
        ax1.tick_params(direction='in', which='both', length=6, width=1, top=True)
        ax1.tick_params(axis='x', labelsize=0)
        ax1.tick_params(axis='y', labelsize=18)
        ax1.set_ylim(bottom=0.1*np.nanmin(min_flux), top=2.*np.nanmax(max_flux))
        ax1.set_xlim(np.nanmin(wave)/1.2, 1.2*np.nanmax(wave))
        ax1.set_ylabel(r'$f_\nu$ (Jy)', fontsize=18)
        #ax1.set_xscale('log')
        ax1.set_yscale('log')

        #xlabs = [1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50, 100, 200, 500]
        #ax1.set_xticks(xlabs[(np.where(xlabs > np.nanmin(wave))[0][0]):(np.where(xlabs < np.nanmax(wave))[0][-1]+1)])
        #ax1.xaxis.set_major_formatter(ScalarFormatter())

        interpMod = np.interp(wave, comps['wave'], fMod)
        res = (flux-interpMod) / flux * 100 # in percentage
        std = np.nanstd(res)
        ax2.plot(wave, res, color='white', linewidth=1)
        ax2.axhline(0., color='white', linestyle='--')
        ax2.tick_params(direction='in', which='both', length=6, width=1,  right=True, top=True)
        ax2.tick_params(axis='x', labelsize=18)
        ax2.tick_params(axis='y', labelsize=14)
        ax2.set_ylim(-4*std, 4*std)
        ax2.set_xlabel(r'$\lambda_{\rm{rest}}$ $(\mu \rm{m})$', fontsize=18)
        ax2.set_ylabel('Residuals (%)', fontsize=14)
        ax2.set_xlim(lamb_low, lamb_high)
        #ax1.set_zorder(100)

        #ax1.set_title('CAFE Spectrum Decomposition', fontsize=16)
        ax1.set_title(f"{self.target}", loc="right", fontsize=20)
        plt.subplots_adjust(hspace=0)
        
        #ax1.set_xlim(12.5, 13.1)
        #ax1.set_ylim(0.18, 2.6)
        
        if save_name is False:
            plt.show()
            return (fig, ax1, ax2)
        else:
            plot_path = os.path.join(self.cafe_output_path, f"{self.target}_SingleExt_r{self.asec}as/{self.target}_SingleExt_r{self.asec}as_linefit.pdf")
            #plot_path = "neiii_unrectified.pdf"
            fig.savefig(plot_path, dpi=1000, format='pdf', bbox_inches='tight')
            plt.close()
    
    
    def line_cutouts(self):
        
        plt.style.use('dark_background')
        gauss = self.recall_gauss()
        extPAH = self.recall_extPAH()
        comps = self.recall_comps()
        spec = self.recall_data()
        wavemod = comps['wave']
        fCir = comps['fCIR']
        fCld = comps['fCLD']
        fCoo = comps['fCOO']
        fWrm = comps['fWRM']
        fHot = comps['fHOT']
        fStb = comps['fSTB']
        fStr = comps['fSTR']
        fDsk = comps['fDSK']
        fLin = comps['fLIN']
        fPAH = comps['fPAH']
        fMod = fCir + fCld + fCoo + fWrm + fHot + fStb + fStr + fDsk + fLin + fPAH
        fCont = fCir + fCld + fCoo + fWrm + fHot + fStb + fStr + fDsk
        
        # Hydrogen lines
        h2_lines = [[], [], [], []]
        h_lines = [[], [], [], []]
        non_h_lines = [[], [], [], []]
        skip_flag = False
        for idx, name in enumerate(gauss[3]):
            if "H2" in name:
                h2_lines[0].append(gauss[0][idx])
                h2_lines[1].append(gauss[1][idx])
                h2_lines[2].append(gauss[2][idx])
                name_split = name.split('_')[0]
                formatted_name = r"H$_2$ 0-0 $S$" + f"({name_split[-1:]})"
                h2_lines[3].append(formatted_name)
            else:
                if "Pfund" in name or "Humphreys" in name:
                    
                    if skip_flag:
                        skip_flag = False
                        continue
                    if idx != len(gauss[3]) - 1:
                        if gauss[3][idx] == gauss[3][idx + 1]:
                            skip_flag = True
                    series = name.split('_')[0][:-2]
                    transition = name.split('_')[0][-2:]
                    formatted_name = series + f" ({transition[0]}-{transition[1]})"
                    h_lines[0].append(gauss[0][idx])
                    h_lines[1].append(gauss[1][idx])
                    h_lines[2].append(gauss[2][idx])
                    h_lines[3].append(formatted_name)
                else:
                    non_h_lines[0].append(gauss[0][idx])
                    non_h_lines[1].append(gauss[1][idx])
                    non_h_lines[2].append(gauss[2][idx])
                    non_h_lines[3].append(name)

        
        fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3)
        plt.subplots_adjust(hspace=0.4)
        plt.subplots_adjust(wspace=0.3)
        fig.set_size_inches(12, 8)
        axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]
        for i in range(9):
            lamb_min = h2_lines[0][i] - 0.025
            lamb_max = h2_lines[0][i] + 0.025
            delta_lamb = lamb_max - lamb_min
            text_x = lamb_min + 0.01 * delta_lamb
            
            min_array = np.absolute(spec['wave'] - lamb_min)
            max_array = np.absolute(spec['wave'] - lamb_max)
            min_index = min_array.argmin()
            max_index = max_array.argmin()
            min2_array = np.absolute(wavemod - lamb_min)
            max2_array = np.absolute(wavemod - lamb_max)
            min2_index = min2_array.argmin()
            max2_index = max2_array.argmin()
            trim_flux = spec['flux'][min_index : max_index]
            trim_fCont = fCont[min2_index : max2_index]
            
            max_y = np.max(trim_flux) * 1.2
            min_y = np.min(trim_fCont) / 1.2
            delta_y = max_y - min_y
            text_y = min_y + 0.90 * delta_y
            
            lflux = gauss_prof(wavemod, [[h2_lines[0][i]], [h2_lines[1][i]], [h2_lines[2][i]]], ext=extPAH)            
            axes[i].plot(wavemod, lflux+fCont, color='yellow', label='_nolegend_', alpha=0.8, linewidth=1, zorder=0)
            #axes[i].plot(wavemod, lflux)
            axes[i].scatter(spec['wave'][min_index : max_index], trim_flux, color="white", s=3, edgecolor='white', facecolor='none', label='Spec Data', alpha=1, zorder=0)
            axes[i].errorbar(spec['wave'][min_index : max_index], spec['flux'][min_index : max_index], yerr=spec['flux_unc'][min_index : max_index], capsize=2, capthick=1, color='white', alpha=1)
            axes[i].plot(wavemod, fCont, color='white', label='Continuum Fit', linestyle="dashed", zorder=2, alpha=0.5)
            axes[i].plot(wavemod, fMod, color='cyan', label='Continuum Fit', linestyle="solid", zorder=2, alpha=0.5)
            axes[i].set_xlim(lamb_min, lamb_max)
            axes[i].set_ylim(min_y, max_y)
            axes[i].text(text_x, text_y, h2_lines[3][i])
            axes[i].set_xlabel(r'$\lambda_{\rm{rest}}$ $(\mu \rm{m})$', fontsize=10)
            axes[i].set_ylabel(r'$f_\nu$ (Jy)', fontsize=10)
            
            fmt = '%.3f'  # as per no of zero(or other) you want after decimal point
            yticks = ticker.FormatStrFormatter(fmt)
            axes[i].yaxis.set_major_formatter(yticks)

        plot_path = os.path.join(self.cafe_output_path, f"{self.target}_SingleExt_r{self.asec}as/{self.target}_h2lines.pdf")
        fig.savefig(plot_path, dpi=1000, format='pdf', bbox_inches='tight')
        plt.close()
        
        
        fig, (ax1, ax2, ax3) = plt.subplots(ncols=3)
        plt.subplots_adjust(wspace=0.3)
        fig.set_size_inches(12, 3)
        axes = [ax1, ax2, ax3]
        
        for i in range(3):
            
            lamb_min = h_lines[0][i] - 0.025
            lamb_max = h_lines[0][i] + 0.025
            delta_lamb = lamb_max - lamb_min
            text_x = lamb_min + 0.01 * delta_lamb
            
            min_array = np.absolute(spec['wave'] - lamb_min)
            max_array = np.absolute(spec['wave'] - lamb_max)
            min_index = min_array.argmin()
            max_index = max_array.argmin()
            min2_array = np.absolute(wavemod - lamb_min)
            max2_array = np.absolute(wavemod - lamb_max)
            min2_index = min2_array.argmin()
            max2_index = max2_array.argmin()
            trim_flux = spec['flux'][min_index : max_index]
            trim_fCont = fCont[min2_index : max2_index]
            
            max_y = np.max(trim_flux) * 1.2
            min_y = np.min(trim_fCont) / 1.2
            delta_y = max_y - min_y
            text_y = min_y + 0.93 * delta_y
            
            lflux = gauss_prof(wavemod, [[h_lines[0][i]], [h_lines[1][i]], [h_lines[2][i]]], ext=extPAH)            
            axes[i].plot(wavemod, lflux+fCont, color='yellow', label='_nolegend_', alpha=0.8, linewidth=1, zorder=0)
            #axes[i].plot(wavemod, lflux)
            axes[i].scatter(spec['wave'][min_index : max_index], trim_flux, color="white", s=3, edgecolor='white', facecolor='none', label='Spec Data', alpha=1, zorder=0)
            axes[i].errorbar(spec['wave'][min_index : max_index], spec['flux'][min_index : max_index], yerr=spec['flux_unc'][min_index : max_index], capsize=2, capthick=1, color='white', alpha=1)
            axes[i].plot(wavemod, fCont, color='white', label='Continuum Fit', linestyle="dashed", zorder=2, alpha=0.5)
            axes[i].plot(wavemod, fMod, color='cyan', label='Continuum Fit', linestyle="solid", zorder=2, alpha=0.5)
            axes[i].set_xlim(lamb_min, lamb_max)
            axes[i].set_ylim(min_y, max_y)
            axes[i].text(text_x, text_y, h_lines[3][i])
            axes[i].set_xlabel(r'$\lambda_{\rm{rest}}$ $(\mu \rm{m})$', fontsize=10)
            axes[i].set_ylabel(r'$f_\nu$ (Jy)', fontsize=10)
            
            fmt = '%.3f'  # as per no of zero(or other) you want after decimal point
            yticks = ticker.FormatStrFormatter(fmt)
            axes[i].yaxis.set_major_formatter(yticks)

        plot_path = os.path.join(self.cafe_output_path, f"{self.target}_SingleExt_r{self.asec}as/{self.target}_hlines.pdf")
        fig.savefig(plot_path, dpi=1000, format='pdf', bbox_inches='tight')
        
        """fig, ax = plt.subplots()
        for i in range(len(gauss[0])):
            if extPAH is None:
                extPAH = np.ones(wavemod.shape)
            lflux = gauss_prof(wavemod, [[gauss[0][i]], [gauss[1][i]], [gauss[2][i]]], ext=extPAH)
            ax.scatter(spec['wave'], spec['flux'], color="white", s=2, edgecolor='black', facecolor='none', label='Spec Data', alpha=1, zorder=0)
            ax.plot(wavemod, lflux, color='#0A31FF', label='_nolegend_', alpha=1, linewidth=1)
            ax.plot(wavemod, fCont)
        plt.show()
        pass"""
    
    
    def relative_velocities(self, wave_c, display=None):
        
        spec_dict = self.recall_data()
        spec_dict["relvel"] = ((const.c * (wave_c - spec_dict["wave"])/spec_dict["wave"]).to(u.kilometer / u.second)).value
        return spec_dict


## TODO
# Second order correction via known line