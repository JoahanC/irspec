""" 
This file defines a datacube object which stores any metadata and real 
data associated with a JWST stage 3 datacube product.
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits 
import astropy.units as u
import astropy.constants as const


class Datacube:
    """ 
    An object which stores metadata and any scientific parameters which
    are associated with a JWST datacube.
    
    Attributes
    ----------
    filename : str
        The location of the JWST datacube.
    general_header : Header
        The non-science high level header of the datacube.
    science_header : Header
        The header containing the observation parameters.
    science_data : ndarray
        The flux array for IFS measurements.
    error_data : ndarray
        The flux error array for IFS measurements.
    dq_data : ndarray
        The data quality flag array for IFS measurements.
    
    Methods
    -------
    
    
    """
    
    
    def __init__(self, filepath, redshift=0, verbose=False):
        """Inits Datacube.
        
        Parameters
        ----------
        filepath : str
            The path to the datacube.
        redshift : float, optional
            The redshift of the observations contained in the datacube.
        verbose : bool, optional
            Whether to display logs.
        """
        
        self.filepath = filepath
        self.redshift = redshift
        self._verbose = verbose
        self._read_data()
        self._read_headers()
        
    
    
    ### Housekeeping Methods
    
    
    def _read_headers(self):
        """Reads and stores header information metadata."""
        self.general_header = fits.getheader(self.filepath, ext=0)
        self.science_header = fits.getheader(self.filepath, ext=1)

        # Housekeeping variables for computations
        self._wv_unit = u.Unit(self.science_header["CUNIT3"])
        self._flux_unit = u.Unit(self.science_header["BUNIT"])
        self._ref_wv = self.science_header["CRVAL3"] * self._wv_unit
        self._delta_wv = self.science_header["CDELT3"] * self._wv_unit
        self.wvs = (self._ref_wv + self._delta_wv * np.arange(self.ax3_len)) / (1 + self.redshift)
        self.area_sr = self.science_header["PIXAR_SR"] * u.steradian
    
    
    def _read_data(self):
        """Reads and stores all relevant data arrays."""
        self.science_data = fits.getdata(self.filepath, ext=1)
        self.error_data = fits.getdata(self.filepath, ext=2)
        self.dq_data = fits.getdata(self.filepath, ext=3)
        self.im_shape = np.shape(self.science_data[0])
        self.ax3_len = len(self.science_data)
    
    
    def display_dq(self, idx):
        """Displays the data quality array for a given array slice."""
        fig, ax = plt.subplots()
        ax.imshow(self.dq_data[idx], origin="lower")
        ax.set_xlabel("XPIX")
        ax.set_ylabel("YPIX")
        plt.show()
    
    
    def spaxel_values(self, x_pix, y_pix):
        """Returns the spaxel flux spectrum, spectrum errors, and 
        spectrum data quality flag at a given (x, y) coordinate."""
        return ((self.science_data[:, y_pix, x_pix] * self._flux_unit * self.area_sr).to(u.Jy),
                (self.error_data[:, y_pix, x_pix] * self._flux_unit * self.area_sr).to(u.Jy),
                self.dq_data[:, y_pix, x_pix])
    
    
    def wv_to_idx(self, wv, unit="um"):
        """Returns the closest index corresponding to a given wavelength
        value.
        
        Parameters
        ----------
        wv : float, Quantity
        unit : Unit
        """
        if not isinstance(wv, u.Quantity):
            wv *= u.Unit(unit)
        wv = wv.to(self.wvs[0].unit)
        if wv >= self.wvs[0] and wv < self.wvs[-1]:
            return np.argmin(abs(self.wvs - wv))
        else:
            raise ValueError("Wavelength is out of bounds!")
    
    
    def wv_to_vel(self, wv, ref_wv, unit="um"):
        """Returns the velocity shift corresponding to a wavelength shift.
        
        Parameters
        ----------
        wv : float, Quantity
            The wavelength at which to evaluate the velocity shift
        ref_wv : float, Quantity
            The reference wavelength to evaluate the wavelength shift.
        unit : Unit
            The units of the provided wavelengths
        
        Returns
        -------
        Quantity : The velocity shift in kilometers/second.
        """
        if not isinstance(wv, u.Quantity):
            wv *= u.Unit(unit)
        if not isinstance(ref_wv, u.Quantity):
            ref_wv *= u.Unit(unit)
        return (const.c * (wv - ref_wv)/ref_wv).to(u.kilometer / u.second)
    
    
    def vel_to_wv(self, vel, ref_wv, vel_unit="km/s", wv_unit="um"):
        """Returns the wavelength corresponding to a velocity shift and 
        a reference wavelength.
        
        Parameters
        ----------
        vel : float, Quantity
            The velocity shift in kilometers/second.
        ref_wv : float, Quantity
            The reference wavelength to evaluate the velocity shift.
        vel_unit : Unit
            The units of the velocity shift.
        wv_unit : Unit
            The units of the reference wavelength and of the returned 
            wavelength.
        
        Returns
        -------
        Quantity : The wavelength shift in the specified `wv_unit`.
        """
        if not isinstance(ref_wv, u.Quantity):
            ref_wv *= u.Unit(wv_unit)
        if not isinstance(vel, u.Quantity):
            vel *= u.Unit(vel_unit)
        return (ref_wv * vel / const.c).to(u.micron)
    
    
    def vel_to_sigma(self, ref_wavelength, vel_disp):
        """_summary_

        Args:
            ref_wavelength (_type_): _description_
            vel_disp (_type_): _description_

        Returns:
            _type_: _description_
        """
        return (ref_wavelength * vel_disp * (u.kilometer / u.second) / const.c).decompose()
    
    
    
    
    
    ### Visualization Methods
    
    
    

