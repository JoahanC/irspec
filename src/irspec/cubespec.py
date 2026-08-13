import os
import glob
import csv
from textwrap import dedent
import numbers

import numpy as np
import astropy.io.fits as fits
import astropy.units as u
import astropy.constants as const
from astropy import coordinates
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from irspec.emission_io import read_line_params, read_line_params3
from specutils import Spectrum1D, SpectralRegion
from specutils.manipulation import extract_region
from lmfit.models import GaussianModel
from specutils.analysis import snr_derived
from irspec.spec_helpers import trim_spec, cut_line, fit_continuum, find_nans
from irspec.fitfuncs import * 
from lmfit import Parameters


# CAFE/CRETA (a local editable install, not on PyPI) and asdf are needed
# only by the extraction/fitting/recall methods; they resolve on first use so
# this module imports without them.
class _LazyImport:
    def __init__(self, module_name, attr=None):
        self._module_name = module_name
        self._attr = attr
        self._target = None

    def _resolve(self):
        if self._target is None:
            import importlib
            try:
                module = importlib.import_module(self._module_name)
            except ImportError as exc:
                raise ImportError(
                    "CubeSpec extraction/fitting requires CAFE/CRETA and asdf."
                    " Install CAFE with: pip install -e <path>/CAFE-master"
                ) from exc
            self._target = getattr(module, self._attr) if self._attr else module
        return self._target

    def __getattr__(self, name):
        return getattr(self._resolve(), name)

    def __call__(self, *args, **kwargs):
        return self._resolve()(*args, **kwargs)


creta = _LazyImport("CRETA.creta")
cafe = _LazyImport("CAFE.cafe")
asdf = _LazyImport("asdf")
pah_drude = _LazyImport("CAFE.component_model", "pah_drude")
gauss_prof = _LazyImport("CAFE.component_model", "gauss_prof")
drude_prof = _LazyImport("CAFE.component_model", "drude_prof")
drude_int_fluxes = _LazyImport("CAFE.component_model", "drude_int_fluxes")


class CubeSpec:
    """ 
    A container for performing spectral extraction and fitting 
    operations with CRETA and CAFE.
    """
    
    def __init__(self, input_path, param_path, param_file, creta_output_path,
                 cafe_output_path, redshift, datacube, mode="AGN",
                 name=None, plot_output_path=None, verbose=False):
        """Initialization method for CubeSpec.

        Args:
            input_path (str): The directory to all of the datacubes to be
                loaded.
            param_path (str): The file path to the parameter file.
            param_file (str): The name of the parameter file.
            creta_output_path (str): The directory to output all CRETA files.
            cafe_output_path (str): The directory to output all CAFE files.
            redshift (float): The redshift at which the datacubes are
                evaluated.
            mode (str): The spectral parameter mode to use CAFE in.
            name (str, optional): The run name used by the runner scripts. When
                set, CAFE outputs are resolved from
                ``{cafe_output_path}/{name}_single/{name}_*`` to match the
                ``{NAME}_single`` layout produced by ``run_cafe.py``. When
                ``None``, the legacy ``{target}_SingleExt_r{r_ap}as`` layout is
                used.
            plot_output_path (str, optional): The directory to write rendered
                plots of this run into (created on demand). When ``None``,
                plots are written alongside the CAFE outputs.
            verbose (bool, optional): Whether to output logs on class
                operations.
        """
        # Set and verify essential path variables
        self.input_path = os.path.abspath(input_path) + '/'
        self.param_path = os.path.abspath(param_path) + '/'
        self.param_file = param_file
        self.param_file_path = os.path.join(self.param_path, param_file)
        self.creta_output_path = os.path.abspath(creta_output_path) + '/'
        self.cafe_output_path = os.path.abspath(cafe_output_path) + '/'
        self.plot_output_path = (os.path.abspath(plot_output_path) + '/'
                                 if plot_output_path is not None else None)
        self.datacube = datacube
        self.line_dict = read_line_params3()
        self._verify_paths()
        
        # Set and verify essential files
        self._verify_files()
        
        # Set essential class fields
        self.c = None
        self.redshift = redshift
        self.mode = mode
        self.run_name = name
        self.c_mod = None
        self.cube_model = None  # populated by perform_grid_fit()
        self.verbose = verbose


    def _result_basename(self):
        """Stem shared by this run's CAFE output files.

        With a run name set this is simply the name (matching ``run_cafe.py``);
        otherwise it is the legacy stem CAFE derives from the extraction file
        name, e.g. ``IR23128-S_SingleExt_r10as`` (aperture with the dot removed).
        """
        if self.run_name is not None:
            return self.run_name
        return f"{self.target}_SingleExt_r{str(self.asec).replace('.', '')}as"


    def _result_dir(self):
        """Directory holding this run's CAFE outputs."""
        if self.run_name is not None:
            return os.path.join(self.cafe_output_path, f"{self.run_name}_single")
        return os.path.join(self.cafe_output_path, self._result_basename())


    def _result_path(self, suffix):
        """Full path to a CAFE output file, e.g. ``suffix='_cafefit.asdf'``."""
        return os.path.join(self._result_dir(), self._result_basename() + suffix)


    def _plot_path(self, suffix):
        """Full path to a rendered plot of this run, e.g. ``suffix='_fit.png'``.

        Plots go to ``plot_output_path`` when one was given (the ``plots/cafe``
        tree managed by ``run_cafe.py``), otherwise they land beside the CAFE
        outputs. The directory is created on demand so plotting never depends
        on the caller having made it first.
        """
        directory = self.plot_output_path or self._result_dir()
        os.makedirs(directory, exist_ok=True)
        return os.path.join(directory, self._result_basename() + suffix)


    
    def _verify_paths(self):
        """Checks if all inputted paths are present prior to extraction
        or fitting and either generates such directories or raises an
        error if a directory containing apriori data is missing.
        """
        if not os.path.exists(self.creta_output_path):
            error_string = f"""
            Specified CRETA output directory:
            {self.creta_output_path}
            does not exist!
            """
            raise NotADirectoryError(dedent(error_string))
        if not os.path.exists(self.cafe_output_path):
            error_string = f"""
            Specified CAFE output directory:
            {self.cafe_output_path}
            does not exist!
            """
            raise NotADirectoryError(dedent(error_string))
        if not os.path.exists(self.input_path):
            error_string = f"""
            Specified CRETA input directory:
            {self.input_path}
            does not exist!
            """
            raise NotADirectoryError(dedent(error_string))
        if not os.path.exists(self.param_file_path):
            error_string = f"""
            Specified parameter file directory:
            {self.param_file_path}
            does not exist!
            """
            raise NotADirectoryError(dedent(error_string))
    
    
    def _verify_files(self):
        """Verifies that inputted data files are present for running this
        pipeline. This method assumes that only the Level3-reduced 
        datacubes are present in the input data folder and that the 
        parameter file files a certain naming structure.
        """
        
        # List all data files
        self.data_files = glob.glob(os.path.join(self.input_path,
                                                 "*.fits"))
        # Validate the datacube count before indexing into the list so that a
        # missing/empty input directory raises a clear error rather than an
        # IndexError. Extraneous files (> 12) are tolerated.
        if len(self.data_files) < 12:
            raise FileNotFoundError(
                f"Datacubes missing! Expected at least 12 *.fits files in "
                f"{self.input_path}, found {len(self.data_files)}."
            )
        test_header = fits.getheader(self.data_files[0], ext=0)
        test_target = test_header["TARGPROP"]
        
        # Check parameter files
        with open(self.param_file_path) as file:
            param_lines = file.readlines()
        self.param_dict = {}
        if "single" in self.param_file_path.split("/")[-1]:
            self.param_dict["type"] = "single"
        if "grid" in self.param_file_path.split("/")[-1]:
            self.param_dict["type"] = "grid"
        self.extraction_id = self.param_file_path.split("_")[1]
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
        
        # Cross check all data files
        for filename in self.data_files:
            current_header = fits.getheader(filename, ext=0)
            current_target = current_header["TARGPROP"]
            if current_target != test_target:
                raise ValueError("Datacubes are of different targets!")
        self.ra = test_header["TARG_RA"]
        self.dec = test_header["TARG_DEC"]
        self.target = test_target
        self.header = test_header
    
    
    def _initialize_CRETA(self):
        """Initializes an instance of CRETA for this particular MIRI
        observation.
        """
        self.c = creta.creta(".")
        return True
    
    
    def rewrite_spec_csv(self):
        """ 
        This method rewrites the CRETA outputted csv into one that is
        accepted by Thomas Lai's web app for viewing spectra.
        """
        if self.verbose:
            method_string = f"""
            Reformating spectrum csv file to be Thomas Lai compliant
            """
            print(dedent(method_string))
        
        # CRETA names its outputs with the raw float radius, e.g. r0.3as
        # (see CRETA/creta.py: 'str(user_rs_arcsec[j])'), so use self.asec
        # directly to match the on-disk file.
        original_csv = os.path.join(self.creta_output_path,
                                    f"{self.target}_SingleExt_r{self.asec}as.csv")
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
    
    
    ### Value recovery routines
    
    
    def open_asdf(self):
        """ 
        This routine opens the asdf file corresponding to the last 
        CAFE fitting session.
        """
        asdf_fn = self._result_path("_cafefit.asdf")
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
                                parfile_path=self.param_path, 
                                parfile_name='/' + self.param_file, 
                                data_path=self.input_path, 
                                output_path=self.creta_output_path,
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
        if self.c == None:
            self._initialize_CRETA()
        if self.param_dict["type"] != "grid":
            raise ValueError("Parameter file is not configured for single extraction.")
        print(self.creta_output_path)
        self.c.gridExtraction(data_path=self.input_path, 
                              parfile_path=self.param_path, 
                              parfile_name=self.param_file, 
                              output_path = self.creta_output_path,
                              output_filebase_name=f'{self.target}',
                              ignore_DQ=True)
    
    
    def read_grid_extraction(self, file):
        if self.c == None:
            self._initialize_CRETA()
        return self.c.customFITSReader(file)
    
    
    def perform_fit(self):
        """ 
        This method runs the CAFE fitting tool on this MIRI observation
        given that CRETA outputs are already available for this object.
        """
        
        # Define appropriate directories and pathing structures.
        # Derive the CAFE package directory from the imported module so this
        # works with an editable install (CAFE is not copied into site-packages).
        cafe_site_dir = os.path.dirname(cafe.__file__) + "/"
        cafe_output_path = self.cafe_output_path #os.path.join(self.cafe_output_path, f"cafe_output/{self.target}/{self.extraction_id}/{self.fit_dirname}")
        source_fd = self.creta_output_path
        source_fn = f'{self.target}_SingleExt_r{self.asec}as.fits'
        
        inppar_fn = cafe_site_dir + f"inp_parfiles/inpars_jwst_nirspec-miri_{self.mode}.ini"
        optpar_fn = cafe_site_dir + "opt_parfiles/default_opt.cafe"
        
        # Initialize the spectroscopic fitting class
        s = cafe.specmod(cafe_site_dir)
        
        # Read in the spectrum
        s.read_spec(source_fn, file_dir=source_fd + '/', z=self.redshift)

        # Thoroughly fit the spectrum
        s.fit_spec(inppar_fn, optpar_fn, output_path=cafe_output_path + '/')


    def _grid_cube_path(self):
        """Locate the CRETA grid-extraction *cube* that CAFE's ``cubemod``
        ingests.

        ``perform_grid_extraction`` writes two products next to each other: a
        per-spaxel *table* (``*_GridExt_*.fits``) and the reshaped *cube*
        (``*_GridExt_*_cube.fits``). CAFE reads the latter. The runner scripts
        additionally rebase these onto the run name, so a cube may also be
        present as ``{run_name}_cube.fits``. Both layouts are supported.
        """
        candidates = []
        if self.run_name is not None:
            candidates.append(f"{self.run_name}_cube.fits")
        candidates.append(f"{self.target}_cube.fits")
        for name in candidates:
            path = os.path.join(self.creta_output_path, name)
            if os.path.exists(path):
                return path
        # Fall back to the raw CRETA naming (size token intact).
        matches = sorted(glob.glob(os.path.join(self.creta_output_path,
                                                "*_GridExt_*_cube.fits")))
        if not matches:
            matches = sorted(glob.glob(os.path.join(self.creta_output_path,
                                                    "*_cube.fits")))
        if not matches:
            raise FileNotFoundError(
                f"No grid-extraction cube (*_cube.fits) found in "
                f"{self.creta_output_path}. Run perform_grid_extraction() first."
            )
        if len(matches) > 1 and self.verbose:
            print(f"Multiple grid cubes found; using {matches[0]}")
        return matches[0]


    def perform_grid_fit(self, extract='Flux_st', trim=True, force_all_lines=False,
                         pattern='snr', resume=False, checkpoint_every=10):
        """PROTOTYPE. Fit every spaxel of a CRETA grid extraction with CAFE.

        This is the grid analogue of :meth:`perform_fit`: it drives CAFE's
        ``cubemod`` instead of ``specmod``. It ingests the grid-extraction cube
        written by :meth:`perform_grid_extraction`
        (``*_GridExt_*_cube.fits``) and fits each grid spaxel individually,
        walking the grid in descending signal-to-noise and seeding every fit
        with the parameters of its highest-SNR already-fitted neighbour (the
        initialisation policy is CAFE's ``fit_cube``).

        CAFE writes a parameter cube (``<stem>_parcube.fits``) and the central
        spaxel's ``<stem>_fitpars.ini`` into ``<cafe_output_path>/<stem>/``.
        Per-parameter 2D maps can then be rendered with
        :meth:`make_grid_maps`.

        Args:
            extract (str, optional): Which CRETA flux column to fit:
                ``'Flux_st'`` (band-stitched, default) or ``'Flux_ap'`` (raw
                aperture flux).
            trim (bool, optional): Drop the overlapping wavelengths between
                adjacent MIRI sub-bands before fitting.
            force_all_lines (bool, optional): Fit every line in the line list
                at every spaxel, including undetected ones. Slower, but yields
                spatially complete maps.
            pattern (str, optional): Spaxel visiting order passed to CAFE's
                ``get_fit_sequence``. ``'snr'`` (default) fits in descending
                signal-to-noise, seeding each spaxel from its highest-SNR
                fitted neighbour. (CAFE's own ``fit_cube`` default of ``None``
                is a bug -- it disables sequence generation -- so a valid mode
                is supplied here.)
            resume (bool, optional): Continue an interrupted fit from its
                checkpoint (``<stem>_parcube.ckpt.fits``, written beside the
                eventual parameter cube) instead of refitting the grid from
                the beginning. Spaxels already fitted are reloaded and still
                seed their neighbours, so the resumed fit follows the same
                path through the grid the uninterrupted one would have. A
                checkpoint from a different cube, mode or redshift is refused
                with a warning rather than mixed into this run.
            checkpoint_every (int, optional): Write that checkpoint every this
                many spaxels. The ceiling on what an interrupt costs is the
                work since the last one plus the spaxel in flight. Set to 0 to
                disable checkpointing entirely.

        Returns:
            CAFE.cafe.cubemod: The fitted cube model, also stored on
            ``self.cube_model``. Its ``.parcube`` holds the per-spaxel best-fit
            parameters.
        """
        # Derive the CAFE package directory from the imported module so this
        # works with an editable install (mirrors perform_fit).
        cafe_site_dir = os.path.dirname(cafe.__file__) + "/"
        cube_path = self._grid_cube_path()
        cube_dir, cube_fn = os.path.split(cube_path)

        inppar_fn = cafe_site_dir + f"inp_parfiles/inpars_jwst_nirspec-miri_{self.mode}.ini"
        optpar_fn = cafe_site_dir + "opt_parfiles/default_opt.cafe"

        if self.verbose:
            print(dedent(f"""
            Fitting grid extraction spaxel-by-spaxel with CAFE
              cube:   {cube_path}
              extract: {extract}
              mode:   {self.mode}
              output: {self.cafe_output_path}
            """))

        # Initialize the cube-fitting class and read the grid cube.
        c = cafe.cubemod(cafe_site_dir)
        c.read_cube(cube_fn, file_dir=cube_dir + '/', extract=extract,
                    trim=trim, z=self.redshift)

        # Fit every spaxel (writes <stem>_parcube.fits + <stem>_fitpars.ini).
        c.fit_cube(inppar_fn, optpar_fn,
                   output_path=self.cafe_output_path + '/',
                   force_all_lines=force_all_lines,
                   pattern=pattern,
                   resume=resume,
                   checkpoint_every=checkpoint_every)

        self.cube_model = c
        return c


    def grid_extract_and_fit(self, extract='Flux_st', force_all_lines=False):
        """PROTOTYPE end-to-end: CRETA grid extraction -> per-spaxel CAFE fit.

        Convenience wrapper that runs :meth:`perform_grid_extraction` and then
        :meth:`perform_grid_fit` on the resulting cube. Requires a grid-type
        parameter file.
        """
        self.perform_grid_extraction()
        return self.perform_grid_fit(extract=extract,
                                     force_all_lines=force_all_lines)


    def make_grid_maps(self, parnames):
        """PROTOTYPE. Write per-parameter 2D maps from the last grid fit.

        ``parnames`` may be a single parameter name or a list (e.g.
        ``['VGRAD', 'COO_TMP']``). Each map is written as
        ``<stem>_<parname>_map.fits`` into the fit's output directory.
        """
        if self.cube_model is None:
            raise RuntimeError(
                "No grid fit available; run perform_grid_fit() first."
            )
        if isinstance(parnames, str):
            parnames = [parnames]
        for parname in parnames:
            self.cube_model.make_map(parname, map_dir=self.cube_model.outPath)


    def generate_extraction_grid(self, gridPointsX, gridPointsY, gridstep_dist, user_ra, user_dec):
        """ 
        This method generates a series of points for doing a series of single extractinos.
        """
        
        NX = np.arange(0,gridPointsX)
        NY = np.arange(0,gridPointsY)
        
        gridstep_dist_pix = gridstep_dist / self.datacube.pixel_scale
        if r_ap == -1:
            print('Warning: The spaxel size is not defined. Using as default the distance between spaxels.') 
            r_ap = gridstep_dist / 2
            r_pix = gridstep_dist_pix / 2
        else:    
            r_pix = r_ap / self.datacube.pixel_scale

        c1 = SkyCoord(user_ra, user_dec, unit="deg")  # defaults to 
        Stringc = SkyCoord(params['user_ra'], params['user_dec'], frame='icrs')     
        user_x, user_y, user_z = self.datacube.wcs.world_to_pixel(c1, self.datacube.wvs[0])
        grids_xs = user_x + (NX - (gridPointsX-1)/2) * gridstep_dist_pix
        grids_ys = user_y + (NY - (gridPointsY-1)/2) * gridstep_dist_pix
        
        
        sky_list = []
        pixels_list = []
        #coord_grid = []
        names = []
        sky_ra = []
        sky_dec = []
        for i in range(len(grids_xs)):
            for j in range(len(grids_ys)):
                sky = self.datacube.pixel_scale.wcs.pixel_to_world(grids_xs[i], grids_ys[j], 0)
                #coord_grid.append(sky)   
                sky_list.append(sky[0])
                pixels_list.append([i,j])
                names.append(str(i)+"_"+str(j))
                sky_ra.append(sky[0].ra)
                sky_dec.append(sky[0].dec)
                
        # for i in range(len(subchannels)):                
        #     self.plotGridSubchanel( user_ra, user_dec, gridstep_dist, gridPointsX, gridPointsY, subchannels[i], r)
        # params_path = current_path+"/Params"
        # self.delteFilesatPath(params_path)
        # self.writeParamsFiles(coord_grid,r,l_ap,pointSource)        
        
        return sky_list, pixels_list, names, sky_ra, sky_dec
    
    
    def return_lines(self):
        return self._result_path("_linetable.ecsv")

    def return_pah(self):
        return self._result_path("_pahtable.ecsv")
    
    def mask_spaxel(self, spaxel, line):
        """Develops the masks necessary to perform gaussian fitting."""
        lower_continuum_region = SpectralRegion(self.line_dict[line][2][0] * self.datacube._wv_unit,
                                                 self.line_dict[line][3][0] * self.datacube._wv_unit)
        upper_continuum_region = SpectralRegion(self.line_dict[line][2][1] * self.datacube._wv_unit,
                                                self.line_dict[line][3][1] * self.datacube._wv_unit)
        continuum_region = lower_continuum_region + upper_continuum_region
        line_region = SpectralRegion(self.line_dict[line][3][0] * self.datacube._wv_unit,
                                        self.line_dict[line][3][1] * self.datacube._wv_unit)
        spaxel_continuum = extract_region(spaxel, continuum_region, return_single_spectrum=True)
        spaxel_line = extract_region(spaxel, line_region)
        return spaxel_continuum, spaxel_line
    
    
    def background_subtraction(self, spaxel_continuum, spaxel_line):
        """Performs background subtraction on a spaxel."""
        
        nanidx = np.where(np.isfinite(spaxel_continuum.flux.value))[0]
        popt, pcov = fit_continuum(spaxel_continuum.spectral_axis.value[nanidx], 
                                   spaxel_continuum.flux.value[nanidx])
        spaxel_continuum_fit = Spectrum1D(OneDPolynomial(spaxel_line.spectral_axis.value, popt[0], popt[1]) * u.Jy, 
                                            spaxel_line.spectral_axis)
        spaxel_continuum_sub = spaxel_line - spaxel_continuum_fit
        line_fluxes = spaxel_continuum_sub.flux.value
        line_wavelengths = spaxel_continuum_sub.spectral_axis.value
        return line_wavelengths, line_fluxes, spaxel_continuum_fit.flux.value
    
    def reconstruct_gaussian_fluxes(self, best_result, spaxel_line, background_flux, line_wavelengths):
        """Regenerates the continuum flux, the total flux, and all of the
        gaussian component fluxes for an lmfit gaussian decomposition model."""
        
        prefixes = ["g1_", "g2_", "g3_"]
        model_components = []
        total_flux = np.copy(background_flux)
        for prefix in prefixes:
            if prefix+"amplitude" in best_result.best_values:
                params = Parameters()
                params.add(name="amplitude", 
                           value=best_result.best_values[prefix + "amplitude"])
                params.add(name="center", 
                           value=best_result.best_values[prefix + "center"])
                params.add(name="sigma", 
                           value=best_result.best_values[prefix + "sigma"])
                model_component = GaussianModel()
                component_flux = model_component.eval(params=params, x=line_wavelengths)
                total_flux += component_flux
                model_components.append(background_flux + component_flux)
        return total_flux, model_components
    
    
    def refit_lines(self, lines):
        """This fits a set of emission lines manually."""
        
        asdf_fn = self._result_path("_cafefit.asdf")
        af = asdf.open(asdf_fn)
        
        wave = np.asarray(af.tree['cafefit']['obsspec']['wave'])
        flux = np.asarray(af['cafefit']['obsspec']['flux'])
        flux_unc = np.asarray(af['cafefit']['obsspec']['flux_unc'])
        
        amplitudes = []
        centers = []
        sigmas = []
        for line in lines:
            nanidx = np.where(np.isfinite(flux))[0]
            spaxel = Spectrum1D(flux[nanidx] * u.Jy, wave[nanidx] * u.micron)
            
            # Process spaxel
            spaxel_continuum, spaxel_line = self.mask_spaxel(spaxel, line)
            line_wavelengths, line_fluxes, background_flux = self.background_subtraction(spaxel_continuum, spaxel_line)
            
            #fig, ax = plt.subplots()
            #ax.plot(spaxel_line.spectral_axis, spaxel_line.flux)
            #ax.plot(spaxel_continuum.spectral_axis, spaxel_continuum.flux)
            #ax.plot(line_wavelengths, line_fluxes)
            #plt.show()
            
            # Estimate line parameters
            line_center = self.line_dict[line][4]
            center_offset = self.datacube.vel_to_wv(self.line_dict[line][5], line_center).value
            center_cutoff = self.datacube.vel_to_wv(self.line_dict[line][6], line_center).value
            minimum_sigma = self.datacube.disp_to_vel(self.line_dict[line][7], line_center).value
            narrow_sigma = self.datacube.disp_to_vel(self.line_dict[line][8], line_center).value
            broad_sigma = self.datacube.disp_to_vel(self.line_dict[line][9], line_center).value
            maximum_sigma = self.datacube.disp_to_vel(self.line_dict[line][10], line_center).value
            narrow_broad_sigma = (minimum_sigma + maximum_sigma) / 2
            minimum_amp = np.sqrt(2*np.pi)*np.max(line_fluxes)*minimum_sigma
            narrow_amp = np.sqrt(2*np.pi)*np.max(line_fluxes)*narrow_sigma
            broad_amp = np.sqrt(2*np.pi)*np.max(line_fluxes)*broad_sigma
            maximum_amp = np.sqrt(2*np.pi)*np.max(line_fluxes)*maximum_sigma * 4
            narrow_broad_amp = (minimum_amp + maximum_amp) / 2
            
            ### Run lmfit fitting routine ###
            
            # Single Gaussian
            single_g1 = GaussianModel(prefix="g1_")
            single_params = single_g1.guess(line_fluxes, x=line_wavelengths)
            single_params.update(single_g1.make_params(center=dict(value=line_center, min=line_center-center_offset, max=line_center+center_offset),
                                        sigma=dict(value=narrow_sigma, min=minimum_sigma, max=maximum_sigma),
                                        amplitude=dict(value=narrow_amp, max=maximum_amp, min=minimum_amp)))
            single_result = single_g1.fit(line_fluxes, single_params, x=line_wavelengths)
            single_redchi = single_result.redchi
            
            # Narrow center + narrow red
            double_g1_2 = GaussianModel(prefix="g1_")
            double_g2_2 = GaussianModel(prefix="g2_")
            double_params_2 = double_g1_2.guess(line_fluxes, x=line_wavelengths)
            double_params_2.update(double_g1_2.make_params(center=dict(value=line_center, min=line_center-center_offset, max=line_center+center_offset),
                                        sigma=dict(value=narrow_sigma, min=minimum_sigma, max=maximum_sigma),
                                        amplitude=dict(value=narrow_amp, max=maximum_amp*10, min=minimum_amp)))
            double_params_2.update(double_g2_2.make_params(center=dict(value=line_center, min=line_center-center_cutoff, max=line_center+center_cutoff),
                                        sigma=dict(value=narrow_sigma, min=minimum_sigma, max=maximum_sigma),
                                        amplitude=dict(value=broad_amp, max=maximum_amp, min=minimum_amp)))
            double_model_2 = double_g1_2 + double_g2_2
            double_result_2 = double_model_2.fit(line_fluxes, double_params_2, x=line_wavelengths)
            
            # Narrow center + narrow blue
            double_g1_3 = GaussianModel(prefix="g1_")
            double_g2_3 = GaussianModel(prefix="g2_")
            double_params_3 = double_g1_3.guess(line_fluxes, x=line_wavelengths)
            double_params_3.update(double_g1_3.make_params(center=dict(value=line_center, min=line_center-center_offset, max=line_center+center_offset),
                                        sigma=dict(value=narrow_sigma, min=minimum_sigma, max=maximum_sigma),
                                        amplitude=dict(value=narrow_amp, max=maximum_amp*10, min=minimum_amp)))
            double_params_3.update(double_g2_3.make_params(center=dict(value=line_center-center_offset, max=line_center, min=line_center-center_cutoff),
                                        sigma=dict(value=narrow_sigma, min=minimum_sigma, max=maximum_sigma),
                                        amplitude=dict(value=broad_amp, max=maximum_amp, min=minimum_amp)))
            double_model_3 = double_g1_3 + double_g2_3
            double_result_3 = double_model_3.fit(line_fluxes, double_params_3, x=line_wavelengths)
            
            # Narrow center + narrow red
            double_g1_4 = GaussianModel(prefix="g1_")
            double_g2_4 = GaussianModel(prefix="g2_")
            double_params_4 = double_g1_4.guess(line_fluxes, x=line_wavelengths)
            double_params_4.update(double_g1_4.make_params(center=dict(value=line_center, min=line_center-center_offset, max=line_center+center_offset),
                                        sigma=dict(value=narrow_sigma, min=minimum_sigma, max=maximum_sigma),
                                        amplitude=dict(value=narrow_amp, max=maximum_amp*10, min=minimum_amp)))
            double_params_4.update(double_g2_4.make_params(center=dict(value=line_center+center_offset, min=line_center, max=line_center+center_cutoff),
                                        sigma=dict(value=broad_sigma, min=minimum_sigma, max=maximum_sigma),
                                        amplitude=dict(value=broad_amp, max=maximum_amp*0.5, min=minimum_amp)))
            double_model_4 = double_g1_4 + double_g2_4
            double_result_4 = double_model_4.fit(line_fluxes, double_params_4, x=line_wavelengths)
            
            # Narrow center + narrow blue
            double_g1_5 = GaussianModel(prefix="g1_")
            double_g2_5 = GaussianModel(prefix="g2_")
            double_params_5 = double_g1_5.guess(line_fluxes, x=line_wavelengths)
            double_params_5.update(double_g1_5.make_params(center=dict(value=line_center, min=line_center-center_offset, max=line_center+center_offset),
                                        sigma=dict(value=narrow_sigma, min=minimum_sigma, max=maximum_sigma),
                                        amplitude=dict(value=narrow_amp, max=maximum_amp*10, min=minimum_amp)))
            double_params_5.update(double_g2_5.make_params(center=dict(value=line_center-center_offset, max=line_center, min=line_center-center_cutoff),
                                        sigma=dict(value=broad_sigma, min=minimum_sigma, max=maximum_sigma),
                                        amplitude=dict(value=broad_amp, max=maximum_amp*0.5, min=minimum_amp)))
            double_model_5 = double_g1_5 + double_g2_5
            double_result_5 = double_model_5.fit(line_fluxes, double_params_5, x=line_wavelengths)
            

            # Record reduced chi squared values
            try:
                single_redchi = single_result.redchi
            except:
                single_redchi = 1e7
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

            if single_redchi <= double_redchi:
                best_fit_type = "single"
            else:
                best_fit_type = "double"
            
            if best_fit_type == "single":
                best_result = single_result
                best_redchi = single_redchi    
            if best_fit_type == "double":
                best_result = double_result
                best_redchi = double_redchi  
            amplitude, center, sigma = self.validate_snr(best_result, 
                                                            spaxel_line, 
                                                            background_flux, 
                                                            line_wavelengths,
                                                            line)
            if len(amplitude) == 0:
                amplitudes.append(0)
            if len(amplitude) == 1:
                amplitudes.append(amplitude[0] * u.watt / u.meter ** 2)
            if len(amplitude) == 2:
                amplitudes.append((amplitude[0] + amplitude[1])* u.watt / u.meter ** 2)
            if len(center) == 0:
                centers.append(0)
            if len(sigma) == 0:
                sigmas.append(0)
        return amplitudes, centers, sigmas
    
    def validate_snr(self, best_result, spaxel_line, background_flux, line_wavelengths, line):
        """Ensures that fitted components pass the SNR threshold and returns
        the array of fitted values."""
        total_flux, model_components = self.reconstruct_gaussian_fluxes(best_result, 
                                                                        spaxel_line, 
                                                                        background_flux, 
                                                                        line_wavelengths
                                                                        )
        amplitudes = []
        centers = []
        sigmas = []
        for idx, component in enumerate(model_components):
            fwhm = best_result.best_values[f"g{idx+1}_sigma"] * 2.3548
            component_center = best_result.best_values[f"g{idx+1}_center"]
            line_region = SpectralRegion((component_center - fwhm) * u.Unit("um"),
                                         (component_center + fwhm) * u.Unit("um"))
            component_spectrum = Spectrum1D(component * u.Unit("Jy"), spaxel_line.spectral_axis)
            amp = best_result.best_values[f"g{idx+1}_amplitude"]
            cen = best_result.best_values[f"g{idx+1}_center"]
            sigma = best_result.best_values[f"g{idx+1}_sigma"]
            rel_vel = self.datacube.wv_to_vel(cen, self.line_dict[line][4])
            vel_disp = ((const.c * (np.abs(sigma))/self.line_dict[line][4]).to(u.kilometer / u.second)).value
            if idx == 0:
                if float(self.line_dict[line][5]) - np.abs(rel_vel.value) < 2:
                    continue
                if float(self.line_dict[line][10]) - np.abs(vel_disp) < 2:
                    continue
            if idx == 1:
                if float(self.line_dict[line][6]) - np.abs(rel_vel.value) < 2:
                    continue
                if float(self.line_dict[line][10]) - np.abs(vel_disp) < 2:
                    continue
            
            if snr_derived(component_spectrum, line_region) > 0.01:
                gamma = np.abs(sigma) * 2.355 / cen
                flux = (1 / np.sqrt(np.pi * np.log(2)) * (np.pi * const.c.to('micron/s') / 2) * (amp * u.Jy * gamma / (cen * u.micron))).to(u.watt / u.meter ** 2)
                rel_vel = self.datacube.wv_to_vel(cen, self.line_dict[line][4])
                vel_disp = ((const.c * (np.abs(sigma))/self.line_dict[line][4]).to(u.kilometer / u.second))#.value
                #amplitudes.append(best_result.best_values[f"g{idx+1}_amplitude"])
                amplitudes.append(flux)
                centers.append(rel_vel)
                sigmas.append(vel_disp)
        return amplitudes, centers, sigmas
    
    def _load_fit(self):
        """Reads this run's fit products out of its ``_cafefit.asdf``.

        Returns:
            tuple: ``(spec_dict, comps, gauss, drude, pahext)``, in the form
            :meth:`cafeplot` expects.

        Raises:
            FileNotFoundError: If this run has not been fitted yet.
        """
        asdf_fn = self._result_path("_cafefit.asdf")
        if not os.path.exists(asdf_fn):
            raise FileNotFoundError(dedent(f"""
            No CAFE fit found for this run:
            {asdf_fn}
            Fit the extraction before plotting it.
            """))
        # ASDF arrays are lazy views onto the open file, so everything is
        # materialised here; the caller outlives the file handle.
        with asdf.open(asdf_fn) as af:
            fit = af.tree['cafefit']
            spec_dict = {key: np.asarray(fit['obsspec'][key])
                         for key in ('wave', 'flux', 'flux_unc')}
            comps = {key: np.asarray(value)
                     for key, value in fit['CompFluxes'].items()}
            extPAH = np.asarray(fit['extComps']['extPAH'])
            gauss = [np.asarray(fit['gauss'][key])
                     for key in ('wave', 'gamma', 'peak')]
            drude = [np.asarray(fit['drude'][key])
                     for key in ('wave', 'gamma', 'peak')]
        return spec_dict, comps, gauss, drude, extPAH


    def plot_fit(self, save=True, suffix="_fit.png", dpi=250, plot_drude=True,
                 title=None, xlim=None, ylim=None, residual='percent'):
        """Renders this run's fitted spectrum from its saved CAFE products.

        Everything plotted is read back from ``{run}_cafefit.asdf``, so an
        existing fit can be re-rendered as many times as needed without
        refitting anything.

        Args:
            save (bool or str): ``True`` writes to ``{plot_dir}/{run}{suffix}``
                (see :meth:`_plot_path`), a string writes to that exact path,
                and ``False`` shows the figure interactively instead.
            suffix (str): File-name suffix used when ``save`` is ``True``.
            dpi (int): Resolution of the saved figure.
            plot_drude (bool): If ``True``, draw the individual PAH Drude
                profiles rather than their sum.
            title (str or None): Title annotating the spectrum panel.
            xlim (tuple or None): ``(left, right)`` wavelength limits in
                microns; defaults to the extent of the observed spectrum.
            ylim (tuple or None): ``(bottom, top)`` flux limits in Jy;
                defaults to a margin around the observed spectrum.
            residual (str): Residual panel units, ``'percent'`` (default) or
                ``'sigma'`` (see :meth:`cafeplot`).

        Returns:
            str or tuple: The path written when saving, otherwise the
            ``(fig, ax1, ax2)`` tuple from :meth:`cafeplot`.
        """
        spec_dict, comps, gauss, drude, extPAH = self._load_fit()

        save_name = self._plot_path(suffix) if save is True else save
        # Assuming there is no phot input.
        # TODO: include phot_dict as input.
        result = self.cafeplot(spec_dict, None, comps, gauss, drude,
                               plot_drude=plot_drude, pahext=extPAH,
                               save_name=save_name, dpi=dpi, title=title,
                               xlim=xlim, ylim=ylim, residual=residual)
        if save_name is False:
            return result
        if self.verbose:
            print(f"Wrote fit figure: {save_name}")
        return save_name


    def plot_channels(self, save=True, suffix="_fit_{slug}.png", dpi=250,
                      plot_drude=True, title=None, channels=None,
                      min_points=10, residual='percent'):
        """Renders the fit as a separate figure per instrument channel.

        Companion to :meth:`plot_fit`: the same fit and the same residual
        panel, but one figure per major channel (NIRSpec and MRS CH1-4, see
        :attr:`CHANNELS`), each covering only that channel's wavelengths, so
        features a few pixels wide in the full-range figure are resolvable.
        Every figure is framed on its own data, which leaves in view the
        components that matter *in that channel*, and shades the +/-1 sigma
        envelope of the spectrum, which at this zoom separates real structure
        from noise.

        Channel ranges are observed-frame; they are de-redshifted here because
        CAFE fits, and this plots, in the rest frame. Channels the extraction
        does not cover are skipped.

        Args:
            save (bool): ``True`` writes one file per channel to
                ``{plot_dir}/{run}{suffix}``; ``False`` shows the figures.
            suffix (str): File-name template for each channel. ``{slug}`` is
                replaced with the channel's short name (``ch1``, ...).
            dpi (int): Resolution of the saved figures.
            plot_drude (bool): If ``True``, draw the individual PAH Drude
                profiles rather than their sum.
            title (str or None): Run description; each figure adds its own
                channel line beneath it.
            channels (sequence or None): ``(label, slug, obs_low, obs_high)``
                entries to render, defaulting to :attr:`CHANNELS`.
            min_points (int): Least spectral points a channel needs before it
                gets a figure.
            residual (str): Residual panel units, ``'percent'`` (default) or
                ``'sigma'`` (see :meth:`cafeplot`).

        Returns:
            list: One entry per rendered channel -- the path written when
            saving, otherwise the ``(fig, ax1, ax2)`` tuple.

        Raises:
            FileNotFoundError: If this run has not been fitted yet.
            ValueError: If the spectrum covers none of the channels.
        """
        spec_dict, comps, gauss, drude, extPAH = self._load_fit()
        wave, flux = spec_dict['wave'], spec_dict['flux']

        results = []
        for label, slug, obs_low, obs_high in (channels or self.CHANNELS):
            low, high = obs_low/(1+self.redshift), obs_high/(1+self.redshift)
            inside = (wave >= low) & (wave <= high)
            if np.count_nonzero(inside) < min_points:
                continue

            # Clip to the data actually present: the bluest and reddest
            # channels usually stop short of their nominal edges.
            xlim = (max(low, np.nanmin(wave[inside])),
                    min(high, np.nanmax(wave[inside])))
            channel_line = (f"{label} ({obs_low:.2f}-{obs_high:.2f} "
                            + r"$\mu$m observed)")
            save_name = self._plot_path(suffix.format(slug=slug)) if save else False
            result = self.cafeplot(
                spec_dict, None, comps, gauss, drude, plot_drude=plot_drude,
                pahext=extPAH, save_name=save_name, dpi=dpi,
                title=channel_line if title is None else f"{title}\n{channel_line}",
                xlim=xlim, ylim=self._log_limits(flux[inside]),
                xscale='linear', legend_loc='best', unc_band=True,
                residual=residual)
            results.append(save_name if save else result)

        if not results:
            raise ValueError(
                "The fitted spectrum covers none of the requested channels "
                f"({', '.join(c[0] for c in (channels or self.CHANNELS))})."
            )
        if save and self.verbose:
            print(f"Wrote {len(results)} channel figures")
        return results


    def recall_fit(self, save=False, **kwargs):
        """
        Recalls the saved fit parameters for a previous CAFE fitting
        session and plots them. Thin wrapper kept for existing callers;
        :meth:`plot_fit` is the full interface.
        """
        return self.plot_fit(save=save, **kwargs)

    
    def line_profiles(self, lines):
        asdf_fn = self._result_path("_cafefit.asdf")
        af = asdf.open(asdf_fn)
        
        wave = np.asarray(af.tree['cafefit']['obsspec']['wave'])
        flux = np.asarray(af['cafefit']['obsspec']['flux'])
        flux_unc = np.asarray(af['cafefit']['obsspec']['flux_unc'])
        
        fig, ax = plt.subplots()
        for line in lines:
            print(self.line_dict[line])
            line_center = self.line_dict[line][4]
            delta_wv = self.datacube.vel_to_wv(1500 * u.kilometer/u.second, line_center)
            low_idx = np.argmin(abs(wave - (line_center - delta_wv.value)))
            high_idx = np.argmin(abs(wave - (line_center + delta_wv.value)))
            cutout_wvs = wave[low_idx:high_idx]
            cutout_vels = []
            for wv in cutout_wvs:
                cutout_vels.append(self.datacube.wv_to_vel(wv, line_center).value)
            #nanidx = np.where(np.isfinite(flux))[0]
            #spaxel = Spectrum1D(flux[nanidx], self.datacube.wvs[nanidx])
            #spaxel_continuum, spaxel_line = self._mask_spaxel(spaxel)
            #line_wavelengths, line_fluxes, background_flux = self.background_subtraction(spaxel_continuum, spaxel_line)
            
            bkg_sub = flux[low_idx:high_idx] - np.min(flux[low_idx:high_idx])
            if len(self.line_dict[line]):
                ax.plot(cutout_vels, bkg_sub / np.max(bkg_sub), label=f"{line}:  {self.line_dict[line][-1]} eV")
            else:
                ax.plot(cutout_vels, bkg_sub / np.max(bkg_sub), label=line)
        ax.legend()
        plt.show()
    
    
    def _mask_spaxel(self, spaxel):
        """Develops the masks necessary to perform gaussian fitting."""
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
    
    def recall_line(self, low_lamb=None, high_lamb=None):
        """ 
        Recalls the saved fit parameters for a previous CAFE fitting 
        session.
        """
        
        asdf_fn = self._result_path("_cafefit.asdf")
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
    
    
    # Continuum components in plotting order: (CompFluxes key, label, colour,
    # line style). Shared by the full-range figure and the per-channel panels
    # so a component looks the same wherever it is drawn.
    _CONT_STYLES = (
        ('fCIR', 'Cirrus',    'tab:cyan',   'solid'),
        ('fCLD', 'Cold',      'tab:blue',   'dashed'),
        ('fCOO', 'Cool',      'green',      'dotted'),
        ('fWRM', 'Warm',      'tab:orange', 'dashdot'),
        ('fHOT', 'Hot',       'brown',      'dashed'),
        ('fSTB', 'Starburst', '#FFEC00',    'dotted'),
        ('fSTR', 'Stellar',   '#FF4500',    'solid'),
        ('fDSK', 'AGN',       'tab:red',    'solid'),
    )

    # Components whose fitted temperature is worth annotating in the legend.
    _TEMP_LABELLED = ('fCLD', 'fCOO', 'fWRM', 'fHOT')

    # Observed-frame coverage of each major channel, read off the input cube
    # headers (CRVAL3 .. CRVAL3 + CDELT3*(NAXIS3-1)), as (label, file slug,
    # low, high). NIRSpec and the four MRS channels each get their own figure
    # from :meth:`plot_channels`.
    CHANNELS = (
        ('NIRSpec G395H', 'nirspec',  2.870,  5.270),
        ('MRS CH1',       'ch1',      4.900,  7.650),
        ('MRS CH2',       'ch2',      7.511, 11.699),
        ('MRS CH3',       'ch3',     11.551, 17.979),
        ('MRS CH4',       'ch4',     17.703, 28.699),
    )


    @staticmethod
    def _model_flux(comps):
        """Total fitted flux: every continuum component plus lines and PAHs."""
        return sum(comps[key] for key, _, _, _ in CubeSpec._CONT_STYLES) \
            + comps['fLIN'] + comps['fPAH']


    @staticmethod
    def _continuum_flux(comps):
        """Summed continuum, i.e. the model without lines and PAHs."""
        return sum(comps[key] for key, _, _, _ in CubeSpec._CONT_STYLES)


    @staticmethod
    def _log_limits(values, pad=0.04):
        """Log-space limits framing ``values`` with a margin of ``pad`` times
        their dynamic range in dex. Returns ``(None, None)`` if nothing is
        positive and finite (so the caller can leave the axis to matplotlib).
        """
        values = np.asarray(values)
        positive = values[np.isfinite(values) & (values > 0)]
        if not positive.size:
            return (None, None)
        lo, hi = np.log10(np.min(positive)), np.log10(np.max(positive))
        margin = pad * max(hi - lo, 1.)
        return (10**(lo - margin), 10**(hi + margin))


    @staticmethod
    def _in_frame(ax, handle):
        """Whether a drawn curve has any point inside ``ax``'s current limits.

        Used to keep legends honest: with the axes framed on the data, a
        component can sit entirely outside the panel, and labelling a curve
        nobody can see is worse than not labelling it.
        """
        if not hasattr(handle, 'get_ydata'):
            return True  # scatter/collection: the data itself, always shown
        x = np.asarray(handle.get_xdata(), dtype=float)
        y = np.asarray(handle.get_ydata(), dtype=float)
        if not x.size or not y.size:
            return True
        (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
        return bool(np.any((x >= x0) & (x <= x1) & (y >= y0) & (y <= y1)))


    @staticmethod
    def _on_continuum(feature, cont, rel_peak=0.05, rel_cont=0.01):
        """A feature profile drawn on the continuum, blanked where it adds
        nothing to it.

        Drude and line profiles are evaluated over the whole wavelength grid,
        so away from their own feature they sit exactly on the continuum.
        Drawn in full, dozens of them stack into an opaque band that hides the
        continuum they are supposed to sit on -- badly so in a zoomed panel.
        Points contributing less than ``rel_peak`` of the profile's own peak
        and less than ``rel_cont`` of the continuum become NaN, leaving each
        profile visible only where it actually shapes the model.
        """
        threshold = np.maximum(rel_peak * np.nanmax(feature), rel_cont * cont)
        return np.where(feature > threshold, feature + cont, np.nan)


    def _draw_fit(self, ax, spec, comps, drude, pahext=None, plot_drude=True,
                  params=None, lw=1., scatter_size=4, unc_band=False):
        """Draws the observed spectrum and every fitted component onto ``ax``.

        Both :meth:`cafeplot` and :meth:`plot_channels` render through this, so
        the two figures cannot drift apart in colour or line style. Only the
        weights scale, via ``lw`` and ``scatter_size``.

        Args:
            ax (matplotlib.axes.Axes): The axis to draw on.
            spec (dict): The observed spectrum (``wave``, ``flux``).
            comps (dict): The component fluxes keyed by component name.
            drude (list): The output parameters of the Drude (PAH) profiles.
            pahext (ndarray or None): The extinction curve applied to the PAHs.
            plot_drude (bool): If ``True``, draw the individual Drude profiles;
                otherwise draw the total PAH contribution.
            params (lmfit.Parameters or None): Fitted parameters, used to
                annotate component temperatures.
            lw (float): Scale factor on every line width.
            scatter_size (float): Marker size of the observed spectrum.
            unc_band (bool): If ``True``, shade the +/-1 sigma envelope of the
                observed spectrum. Worth having in a zoomed panel, where it
                separates real structure from noise; over the full range the
                band is narrower than the data trace and only thickens it.

        Returns:
            tuple: ``(wavemod, fCont, fMod)``.
        """
        wavemod = comps['wave']
        fCont = self._continuum_flux(comps)
        fMod = fCont + comps['fLIN'] + comps['fPAH']
        alpha = 1

        unc = spec.get('flux_unc') if hasattr(spec, 'get') else None
        if unc_band and unc is not None and np.any(np.isfinite(unc)):
            ax.fill_between(spec['wave'], spec['flux']-unc, spec['flux']+unc,
                            color='white', alpha=0.25, linewidth=0, zorder=-1,
                            label=r'Spec $\pm1\sigma$')

        ax.scatter(spec['wave'], spec['flux'], color="black", s=scatter_size,
                   linewidths=0.5*lw, edgecolor='white', facecolor='black',
                   label='Spec Data', alpha=0.6, zorder=0)
        ax.plot(wavemod, fCont, color='red', label='Continuum Fit',
                linestyle="dashed", zorder=4, lw=1.2*lw, alpha=0.8)
        ax.plot(wavemod, fMod, color='cyan', label='Total Fit',
                linewidth=lw, zorder=0, alpha=alpha)

        for key, label, color, style in self._CONT_STYLES:
            if not np.any(comps[key] > 0):
                continue
            if params is not None and key in self._TEMP_LABELLED:
                label += r' ({:.0f}$\,$K)'.format(params[key[1:]+'_TMP'].value)
            ax.plot(wavemod, comps[key], label=label, c=color, alpha=alpha,
                    linewidth=lw, linestyle=style)

        if np.any(comps['fLIN'] > 0):
            ax.plot(wavemod, self._on_continuum(comps['fLIN'], fCont, rel_peak=0),
                    label='Lines', c='green', alpha=alpha, linewidth=lw, zorder=1)

        # PAH features, either resolved into their Drude profiles or summed.
        if plot_drude is True:
            if pahext is None:
                pahext = np.ones(wavemod.shape)
            for i in range(len(drude[0])):
                dflux = drude_prof(wavemod, [[drude[0][i]], [drude[1][i]], [drude[2][i]]], ext=pahext)
                # Every profile carries the label, not just the first: the
                # legend keeps one key per label, and in a zoomed channel the
                # first Drude is often not the one on screen.
                ax.plot(wavemod, self._on_continuum(dflux, fCont), color='fuchsia',
                        label='PAHs', alpha=alpha, linewidth=0.7*lw)
        elif np.any(comps['fPAH'] > 0):
            ax.plot(wavemod, fCont+comps['fPAH'], label='PAHs', color='fuchsia',
                    alpha=alpha, linewidth=lw)

        return wavemod, fCont, fMod


    def cafeplot(self, spec, phot, comps, gauss, drude, vgrad={'VGRAD':0.}, plot_drude=True, pahext=None, save_name=False, params=None, dpi=250, title=None, xlim=None, ylim=None, xscale='log', legend_loc='upper left', residual='percent', residual_ylim=None, unc_band=False):
        """Plots the SED and the CAFE fit.

        Args:
            spec (dict): The observed spectrum, with ``wave`` (rest
                wavelength), ``flux``, and ``flux_unc`` entries.
            phot (dict or None): The photometry to overplot, or ``None`` to
                plot the spectrum alone.
            comps (dict): The component fluxes keyed by component name.
            gauss (list): The output parameters of the Gaussian line profiles.
            drude (list): The output parameters of the Drude (PAH) profiles.
            vgrad (dict): The velocity-gradient parameter used to shift the
                model (default ``{'VGRAD': 0.}``).
            plot_drude (bool): If ``True``, plot individual Drude profiles;
                otherwise plot the total PAH contribution (default ``True``).
            pahext (ndarray or None): If not ``None``, the extinction curve
                applied to the PAHs.
            save_name (bool or str): If ``False``, show the figure; otherwise
                the file path to save it to.
            params (lmfit.Parameters or None): The fitted parameters, used to
                annotate component temperatures (default ``None``).
            dpi (int): Resolution of the saved figure (default ``250``).
            title (str or None): Title annotating the spectrum panel.
            xlim (tuple or None): ``(left, right)`` wavelength limits in
                microns; defaults to the extent of the observed spectrum.
            ylim (tuple or None): ``(bottom, top)`` flux limits in Jy;
                defaults to a margin around the observed spectrum.
            xscale (str): ``'log'`` for the full range, where it keeps the
                short wavelengths readable, or ``'linear'`` for a zoom into
                one channel, where a decade-based axis has nothing to say.
            legend_loc (str): Legend placement passed to matplotlib.
            residual (str): ``'percent'`` (default) plots the residual as a
                percentage of the flux. ``'sigma'`` divides by the measurement
                error instead and adds +/-1 and +/-3 sigma guides -- only
                meaningful once the quoted uncertainties are trustworthy. On
                these extractions they are not: they imply SNR ~1700 per
                channel while the fit residuals run ~2%, so sigma residuals
                come out in the tens to hundreds and describe model mismatch,
                not noise. A spectrum with no usable uncertainties falls back
                to percent regardless.
            residual_ylim (tuple or None): ``(bottom, top)`` for the residual
                panel. Defaults to +/-19% for percentage residuals, and for
                sigma residuals to the 99th percentile of |residual| in view,
                never tighter than the +/-3 sigma guide.
            unc_band (bool): If ``True``, shade the +/-1 sigma envelope of the
                observed spectrum on the flux panel.

        Returns:
            tuple or None: ``(fig, ax1, ax2)`` when the figure is shown,
            ``None`` when it is saved (the figure is closed).
        """
        
        plt.style.use('dark_background')
        wavemod = comps['wave']
        fMod = self._model_flux(comps)

        fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios':[3,1]}, figsize=(6,15), sharex=True)
        self._draw_fit(ax1, spec, comps, drude, pahext=pahext,
                       plot_drude=plot_drude, params=params, unc_band=unc_band)
        #ax1.errorbar(spec['wave'], spec['flux'], yerr=spec['flux_unc'], fmt='none', color='white', alpha=0.1)
        if phot is not None:
            ax1.scatter(phot['wave'], phot['flux'], marker='x', s=18, edgecolor='none', facecolor='black', label='Phot Data', alpha=0.9)
            ax1.errorbar(phot['wave'], phot['flux'], xerr=phot['width']/2, yerr=phot['flux_unc'], fmt='none', color='white', alpha=0.1)
            wave = np.concatenate((spec['wave'], phot['wave']))
            flux = np.concatenate((spec['flux'], phot['flux']))
            sortinds = np.argsort(wave)
            unc = np.concatenate((spec['flux_unc'], phot['flux_unc']))[sortinds]
            wave = wave[sortinds] ; flux = flux[sortinds]
        else:
            wave = spec['wave']
            flux = spec['flux']
            unc = spec.get('flux_unc') if hasattr(spec, 'get') else None

        #ax11 = ax1.twinx()
        #ax11.plot(wavemod, pahext, linestyle='dashed', color='pink', alpha=1, linewidth=1, label="Attenuation")
        #ax11.set_ylim(0, 1.1)
        #ax11.set_ylabel(r'Attenuation fraction $_{\rm{Warm\,dust, PAHs, Lines}}$', fontsize=14)
        #ax11.tick_params(axis='y', labelsize=10)
        #ax11.tick_params(direction='in', which='both', length=4, width=0.8, right=True)

        ax1.tick_params(direction='in', which='both', length=6, width=1, top=True)
        ax1.tick_params(axis='x', labelsize=0)
        ax1.tick_params(axis='y', labelsize=16)
        if ylim is None:
            # Frame the observed spectrum itself. Continuum components far
            # below the data then fall out of frame, which is the point:
            # anything that leaves was not a meaningful contributor here.
            ylim = self._log_limits(flux)
        ax1.set_ylim(bottom=ylim[0], top=ylim[1])
        if xlim is None:
            xlim = (np.nanmin(wave)/1.2, 1.2*np.nanmax(wave))
        ax1.set_xlim(left=xlim[0], right=xlim[1])
        if title is not None:
            ax1.set_title(title, loc="right", fontsize=16)

        # Legend last, so components that fell out of the (data-driven) frame
        # are dropped from it rather than labelling a curve that isn't there.
        shown, seen = [], set()
        for handle, label in zip(*ax1.get_legend_handles_labels()):
            if label not in seen and self._in_frame(ax1, handle):
                seen.add(label)
                shown.append((handle, label))
        handles, labels = zip(*shown) if shown else ((), ())
        ax1.legend(handles, labels, loc=legend_loc, fontsize=11,
                   framealpha=0.65, labelspacing=0.3)
        ax1.set_ylabel(r'$f_\nu$ (Jy)', fontsize=20)
        ax1.xaxis.grid(True, which='major', alpha=0.25, color="white")
        ax1.yaxis.grid(True, which='major', alpha=0.25, color="white")
        #ax1.set_xscale('log')
        ax1.set_yscale('log')
        #ax1.axvline(9.7, linestyle='--', alpha=0.2)

        interpMod = np.interp(wave, comps['wave'], fMod)
        usable_unc = (unc is not None
                      and np.any(np.isfinite(unc) & (np.asarray(unc) > 0)))
        if residual == 'sigma' and usable_unc:
            # Data minus model in units of the measurement error: unlike a
            # percentage, this is comparable between bright and faint
            # channels, so the panel reads as significance.
            unc = np.asarray(unc, dtype=float)
            with np.errstate(divide='ignore', invalid='ignore'):
                res = np.where(unc > 0, (flux-interpMod)/unc, np.nan)
            ax2.set_ylabel(r'Residuals ($\sigma$)', fontsize=20)
            guides = (1., 3.)
        else:
            res = (flux-interpMod) / flux * 100 # in percentage
            ax2.set_ylabel('Residuals (%)', fontsize=20)
            guides = ()

        ax2.plot(wave, res, color='white', linewidth=0.7, zorder=3)
        ax2.axhline(0., color='black', linestyle='--')
        ax2.tick_params(direction='in', which='both', length=6, width=1,  right=True, top=True)
        ax2.tick_params(axis='x', labelsize=16)
        ax2.tick_params(axis='y', labelsize=16)

        if residual_ylim is not None:
            span = None
            ax2.set_ylim(*residual_ylim)
        elif guides:
            # Scale to the bulk of the residuals rather than their worst
            # spikes, but never crop tighter than the +/-3 sigma guide.
            visible = res[(wave >= xlim[0]) & (wave <= xlim[1])] if xlim else res
            span = 3.5
            if np.any(np.isfinite(visible)):
                span = max(span, float(np.nanpercentile(np.abs(visible), 99)))
            ax2.set_ylim(-span, span)
        else:
            span = 19.
            ax2.set_ylim(-19, 19)
        for guide in guides:
            if span is None or guide < span:
                for sign in (-1, 1):
                    ax2.axhline(sign*guide, color='white', linestyle=':',
                                alpha=0.35, linewidth=0.8)

        ax2.set_xlabel(r'$\lambda_{\rm{rest}}$ $(\mu \rm{m})$', fontsize=20)
        #ax1.set_zorder(100)
        ax1.set_xscale(xscale)

        if xscale == 'log':
            # A log axis otherwise labels only the decade, which leaves an MRS
            # spectrum with a single tick; label the round wavelengths in range
            # instead (set after set_xscale, which resets locators/formatters).
            xlabs = np.array([1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50, 100, 200, 500])
            ticks = xlabs[(xlabs >= xlim[0]) & (xlabs <= xlim[1])]
            if ticks.size:
                ax1.set_xticks(ticks)
                ax1.xaxis.set_major_formatter(ScalarFormatter())
                ax1.xaxis.set_minor_formatter(ticker.NullFormatter())

        fig.set_size_inches(10, 6)
        plt.subplots_adjust(hspace=0)
        # Keep the residual label off the flux panel's lowest tick label.
        fig.align_ylabels([ax1, ax2])

        if save_name is False:
            plt.show()
            return (fig, ax1, ax2)
        else:
            os.makedirs(os.path.dirname(os.path.abspath(save_name)), exist_ok=True)
            fig.savefig(save_name, dpi=dpi, bbox_inches='tight')
            plt.close(fig)
            return None
    
    
    def line_diagnostic(self, spec, phot, comps, gauss, drude, line_names, line_waves, lamb_low=None, lamb_high=None, vgrad={'VGRAD':0.}, plot_drude=True, pahext=None, save_name=False, params=None ):
        """Plots the SED and all fitted emission lines for a given target
        within some wavelength range.

        Args:
            lamb_low (float): The lower limit of wavelength to be plotted.
            lamb_high (float): The upper limit of wavelength to be plotted.
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
            plot_path = self._result_path("_linefit.pdf")
            #plot_path = "neiii_unrectified.pdf"
            fig.savefig(plot_path, dpi=1000, format='pdf', bbox_inches='tight')
            plt.close()
    
    
    def line_cutouts(self):
        """ 
        
        """
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

        plot_path = self._result_path("_h2lines.pdf")
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

        plot_path = self._result_path("_hlines.pdf")
        fig.savefig(plot_path, dpi=1000, format='pdf', bbox_inches='tight')
    
    
    def relative_velocities(self, wave_c, display=None):
        
        spec_dict = self.recall_data()
        spec_dict["relvel"] = ((const.c * (wave_c - spec_dict["wave"])/spec_dict["wave"]).to(u.kilometer / u.second)).value
        return spec_dict


## TODO
# Second order correction via known line