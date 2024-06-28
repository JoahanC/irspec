from datacube import Datacube
from pathlib import Path
import astropy.io.fits as fits
from astropy import coordinates
import astropy.units as u
from astroquery.ipac.ned import Ned
import CRETA.creta as creta
import CAFE.cafe as cafe
import numpy as np
import os 
import glob
import site


class CubeSpec:
    """ 
    A container for performing spectral extraction and fitting 
    operations with CRETA and CAFE.
    """
    
    def __init__(self, creta_dir, param_path, param_file, data_path, full_set=True, redshift=None):
        
        self._creta_dir = os.path.abspath(creta_dir)
        self._param_path = os.path.join(self._creta_dir, param_path)
        self.param_path = os.path.join(self._param_path, param_file)
        self.param_file = param_file
        self.data_path = os.path.join(self._creta_dir, data_path)
        self.c = None
        if full_set:
            self._verify_files()
        self.output_path = os.path.join(self._creta_dir, f"creta_output/extractions/{self.target}")
        self._verify_paths()
        self.redshift = redshift
        if self.redshift == None:
            self.find_redshift()
    
    
    def _verify_files(self):
        
        self.data_files = glob.glob(os.path.join(self.data_path, "*.fits"))
        test_header = fits.getheader(self.data_files[0], ext=0)
        test_target = test_header["TARGPROP"]
        if len(self.data_files) == 12:
            for filename in self.data_files:
                current_header = fits.getheader(filename, ext=0)
                current_target = current_header["TARGPROP"]
                if current_target != test_target:
                    raise ValueError("Datacubes are of different targets!")
            self.target = test_target
            self.header = test_header
        else:
            if len(self.data_files) < 12:
                raise FileNotFoundError("Datacubes missing!")
            if len(self.data_files) > 12:
                raise FileExistsError("Extraneous files found!")
    
    
    def _verify_paths(self):
        """ 
        Checks if all inputted paths are present prior to extraction 
        or fitting and either generates such directories or raises an 
        error if a directory containing apriori data is missing.
        """
        if not os.path.exists(self.output_path):
            os.makedirs(self.output_path)
    
    
    def _initialize_CRETA(self):
        """ 
        Verifies the existence of a parameter file and initializes an
        instance of CRETA for this particular MIRI observation.
        """
        pf = open(os.path.join(os.path.join(self._creta_dir, self._param_path), self.param_file), 'r')
        print(pf.read())
        pf.close()
        self.c = creta.creta(self._creta_dir)
    
    
    def perform_extraction(self):
        if self.c == None:
            self._initialize_CRETA()
        self.c.singleExtraction(parameter_file=True, parfile_path=self._param_path, parfile_name="/" + self.param_file, data_path=self.data_path, output_path=self.output_path + "/", output_filebase_name=f'{self.target}')
    
    
    def rewrite_spec_csv(self):
        """ 
        This method rewrites the CRETA outputted csv into one that is 
        accepted by Thomas Lai's web app for viewing spectra.
        """
        import csv 
        
        original_csv = os.path.join(self.output_path, f"{self.target}_SingleExt_r0.7as.csv")
        spec_dict = {"w": [], "f": [], "f_unc": []}
        readlines = False
        with open(original_csv, 'r') as csvfile:
            # Read lines until results appear
            for line in csvfile.readlines():
                print(line)
                if "Wave,Band_name,Flux_ap,Err_ap,R_ap,Flux_ap_st,Err_ap_st,DQ" in line:
                    readlines = True 
                    continue
                if readlines:
                    vals = line.split(sep=",")
                    print(vals)
                    spec_dict["w"].append(vals[0])
                    spec_dict["f"].append(vals[5])
                    spec_dict["f_unc"].append(vals[6])
    
        new_csv = os.path.join(self.output_path, f"{self.target}_sum_sf_spec.csv")
        with open(new_csv, 'w') as csvfile:
            fieldnames = ["w", "f", "f_unc"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for idx, _ in enumerate(spec_dict["w"]):
                writer.writerow({"w": spec_dict["w"][idx], "f" : spec_dict["f"][idx], "f_unc" : spec_dict["f_unc"][idx]})
    
    
    def find_redshift(self):
        print("Autoquerying for redshift using header values.")
        print(f"RA: {self.header['TARG_RA']}, DEC: {self.header['TARG_DEC']}")
        co = coordinates.SkyCoord(ra=self.header["TARG_RA"], dec=self.header["TARG_DEC"], unit=(u.deg, u.deg))
        result_table = Ned.query_region(co, radius=2 * u.arcsec)
        print(f"Found {len(result_table[~np.isnan(result_table['Redshift'])])} result: z = {result_table['Redshift'][0]}")
        self.redshift = result_table["Redshift"]
    
    
    def perform_fit(self):
        cafe_dir = site.getsitepackages()[0]+'/CAFE/'
        creta_dir = site.getsitepackages()[0]+'/CRETA/'
        redshift = 0.055206
        data_path = self.output_path
        cafe_output_path = os.path.join(self._creta_dir, f"cafe_output/{self.target}")
        if not os.path.exists(data_path):
            os.makedirs(data_path)
        else:
            print("Input Data folder already exists.")

        if not os.path.exists(cafe_output_path):
            os.makedirs(cafe_output_path)
        else:
            print("CAFE Output folder already exists.")
        
        source_fd = data_path
        source_fn = f'/{self.target}_SingleExt_r0.7as.fits'

        inppar_fn = cafe_dir+'inp_parfiles/inpars_jwst_miri_AGN.ini'
        optpar_fn = cafe_dir+'opt_parfiles/default_opt.cafe'
        
        s = cafe.specmod(cafe_dir)
        
        # Read in the spec
        s.read_spec(source_fn, file_dir=source_fd, z=self.redshift)

        # Preview the spectrum and overlay the initial params
        s.plot_spec_ini(inppar_fn, optpar_fn)
        
        s.fit_spec(inppar_fn, optpar_fn, output_path=cafe_output_path)
        
        s.plot_spec_fit(inppar_fn, optpar_fn)
        self.cafeplot()
    
    
    def run_spectool(self):
        """ 
        This script runs Thomas Lai's spectool and can be used to make
        bird's eye observations of a given spectra and any stand-out 
        line features.
        """
        import subprocess 
        
        os.chdir("./../JWST-SpecTool/src/")
        
        subprocess.run(["python", "app.py", f"{self.redshift}", f"{self.target}"])
        subprocess.run(["python", "-m", "webbrowser", "-t", "http://127.0.0.1:8050"])
    
    
    
    """def cafeplot(self, spec, phot, comps, gauss, drude, vgrad={'VGRAD':0.}, plot_drude=True, pahext=None, save_name=False, params=None):
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
        ax1.scatter(spec['wave'], spec['flux'], marker='o', s=6, edgecolor='k', facecolor='none', label='Spec Data', alpha=0.9)
        ax1.errorbar(spec['wave'], spec['flux'], yerr=spec['flux_unc'], fmt='none', color='k', alpha=0.1)
        if phot is not None:
            ax1.scatter(phot['wave'], phot['flux'], marker='x', s=18, edgecolor='none', facecolor='k', label='Phot Data', alpha=0.9)
            ax1.errorbar(phot['wave'], phot['flux'], xerr=phot['width']/2, yerr=phot['flux_unc'], fmt='none', color='k', alpha=0.1)
            wave = np.concatenate((spec['wave'], phot['wave']))
            flux = np.concatenate((spec['flux'], phot['flux']))
            sortinds = np.argsort(wave)
            wave = wave[sortinds] ; flux = flux[sortinds]
        else:
            wave = spec['wave']
            flux = spec['flux']
                                
        ax1.plot(wavemod, fCont, color='gray', label='Continuum Fit', linestyle='-', zorder=4, alpha=0.8)
        ax1.plot(wavemod, fCont+fLin+fPAH, color='#4c956c', label='Total Fit', linewidth=1.5, zorder=5, alpha=0.85) # green

        CLD_TMP = '' if params == None else r' ('+"{:.0f}".format(params['CLD_TMP'].value)+'$\,$K'+')'
        COO_TMP = '' if params == None else r' ('+"{:.0f}".format(params['COO_TMP'].value)+'$\,$K'+')'
        WRM_TMP = '' if params == None else r' ('+"{:.0f}".format(params['WRM_TMP'].value)+'$\,$K'+')'
        HOT_TMP = '' if params == None else r' ('+"{:.0f}".format(params['HOT_TMP'].value)+'$\,$K'+')'
            
        alpha = 0.6
        lw = 0.8
        if np.any(fCir > 0):
            ax1.plot(wavemod, fCir, label='Cirrus', c='tab:cyan', alpha=alpha, linewidth=lw)
        if np.sum(fCld > 0):
            ax1.plot(wavemod, fCld, label='Cold'+CLD_TMP, c='tab:blue', alpha=alpha, linewidth=lw)
        if np.any(fCoo > 0):
            ax1.plot(wavemod, fCoo, label='Cool'+COO_TMP, c='#008080', alpha=alpha, linewidth=lw) # teal
        if np.any(fWrm > 0):
            ax1.plot(wavemod, fWrm, label='Warm'+WRM_TMP, c='tab:orange', alpha=alpha, linewidth=lw)
        if np.any(fHot > 0):
            ax1.plot(wavemod, fHot, label='Hot'+HOT_TMP, c='#FFD700', alpha=alpha, linewidth=lw) # gold
        if np.any(fStb > 0): 
            ax1.plot(wavemod, fStb, label='Starburst', c='tab:brown', alpha=alpha, linewidth=lw)
        if np.any(fStr > 0):
            ax1.plot(wavemod, fStr, label='Stellar', c='#FF4500', alpha=alpha, linewidth=lw) # orangered
        if np.any(fDsk > 0):
            ax1.plot(wavemod, fDsk, label='AGN', c='tab:red', alpha=alpha, linewidth=lw)
        if np.any(fLin > 0):
            ax1.plot(wavemod, fCont+fLin, label='Lines', c='#1e6091', alpha=alpha, linewidth=lw) # blue

        # Plot lines
        for i in range(len(gauss[0])):
            if pahext is None:
                pahext = np.ones(wavemod.shape)
            lflux = gauss_prof(wavemod, [[gauss[0][i]], [gauss[1][i]], [gauss[2][i]]], ext=pahext)
            
            ax1.plot(wavemod, lflux+fCont, color='#1e6091', label='_nolegend_', alpha=alpha, linewidth=0.4)
            #if i == 0:
            #    ax1.plot(wavemod, lflux+fCont, color='#1e6091', label='Lines', alpha=alpha, linewidth=0.4)
            #else:
            #    ax1.plot(wavemod, lflux+fCont, color='#1e6091', label='_nolegend_', alpha=alpha, linewidth=0.4)
        
        # Plot PAH features
        if plot_drude is True:
            for i in range(len(drude[0])):
                if pahext is None: 
                    pahext = np.ones(wavemod.shape)
                dflux = drude_prof(wavemod, [[drude[0][i]], [drude[1][i]], [drude[2][i]]], ext=pahext)

                if i == 0:
                    ax1.plot(wavemod, dflux+fCont, color='purple', label='PAHs', alpha=alpha, linewidth=0.5)
                else:
                    ax1.plot(wavemod, dflux+fCont, color='purple', label='_nolegend_', alpha=alpha, linewidth=0.5)
        elif np.any(fPAH > 0):
            ax1.plot(wavemod, fCont+fPAH, label='PAHs', color='purple', alpha=alpha)

        ax11 = ax1.twinx()
        ax11.plot(wavemod, pahext, linestyle='dashed', color='gray', alpha=0.5, linewidth=0.6)
        ax11.set_ylim(0, 1.1)
        ax11.set_ylabel(r'Attenuation fraction $_{\rm{Warm\,dust, PAHs, Lines}}$', fontsize=14)
        ax11.tick_params(axis='y', labelsize=10)
        #ax11.tick_params(direction='in', which='both', length=4, width=0.8, right=True)

        min_flux = np.nanmin(spec['flux'][np.r_[0:5,-5:len(spec['flux'])]])
        max_flux = np.nanmax(spec['flux'][np.r_[0:5,-5:len(spec['flux'])]])

        ax1.legend(loc='lower right')
        ax1.tick_params(direction='in', which='both', length=6, width=1, top=True)
        ax1.tick_params(axis='x', labelsize=0)
        ax1.tick_params(axis='y', labelsize=12)
        ax1.set_ylim(bottom=0.1*np.nanmin(min_flux), top=2.*np.nanmax(max_flux))
        #ax1.set_xlim(left=2.5, right=36)
        ax1.set_xlim(np.nanmin(wave)/1.2, 1.2*np.nanmax(wave))
        ax1.set_ylabel(r'$f_\nu$ (Jy)', fontsize=14)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        #ax1.axvline(9.7, linestyle='--', alpha=0.2)

        xlabs = [1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50, 100, 200, 500]
        ax1.set_xticks(xlabs[(np.where(xlabs > np.nanmin(wave))[0][0]):(np.where(xlabs < np.nanmax(wave))[0][-1]+1)])
        ax1.xaxis.set_major_formatter(ScalarFormatter())

        interpMod = np.interp(wave, comps['wave'], fMod)
        res = (flux-interpMod) / flux * 100 # in percentage
        std = np.nanstd(res)
        ax2.plot(wave, res, color='k', linewidth=1)
        #ax2.plot(wave, (spec['flux']-interpMod)/func, color='k')
        ax2.axhline(0., color='k', linestyle='--')
        ax2.tick_params(direction='in', which='both', length=6, width=1,  right=True, top=True)
        ax2.tick_params(axis='x', labelsize=12)
        ax2.tick_params(axis='y', labelsize=12)
        ax2.set_ylim(-4*std, 4*std)
        #ax2.set_ylim(bottom=-4, top=4)
        ax2.set_xlabel(r'$\lambda_{\rm{rest}}$ $(\mu \rm{m})$', fontsize=14)
        #ax2.set_ylabel(r'$f^{data}_\nu - f^{tot}_\nu$ $(\sigma)$', fontsize=14)
        ax2.set_ylabel('Residuals (%)', fontsize=14)

        #ax1.set_zorder(100)

        #ax1.set_title('CAFE Spectrum Decomposition', fontsize=16)
        plt.subplots_adjust(hspace=0)
        
        # Use black as the patch backgound 
        # fig.patch.set_facecolor('k')
        # ax1.xaxis.label.set_color('w')
        # ax1.yaxis.label.set_color('w')
        # ax1.tick_params(direction='out', which='both', axis='both', colors='w')
        # ax11.tick_params(direction='out', which='both', length=4, width=0.8, right=True, colors='w')
        # ax11.yaxis.label.set_color('w')
        # ax2.xaxis.label.set_color('w')
        # ax2.yaxis.label.set_color('w')
        # #ax11.tick_params(axis='both', colors='w')
        # ax2.tick_params(direction='out', which='both', axis='both', colors='w')
        
        if save_name is False:
            plt.show()
            return (fig, ax1, ax2)
        else:
            fig.savefig(save_name, dpi=500, format='png', bbox_inches='tight')
            plt.close()"""


## TODO
# Second order correction via known line