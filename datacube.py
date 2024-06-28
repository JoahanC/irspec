""" 
This file defines a datacube object which stores any metadata and real 
data associated with a JWST stage 3 datacube product.
"""
import numpy as np
from astropy.io import fits 
import astropy.units as u


class Datacube:
    """ 
    An object which stores metadata and any scientific parameters which
    are associated with a JWST datacube.
    
    Arguments
    ---------
    filename : str
        The location of the JWST datacube.
    """
    
    
    def __init__(self, filename):
        
        
        self.filename = filename
        
        self._read_headers()
        self._read_data()
    
    
    def _read_headers(self):
        """ 
        Runs a routine which reads relevant header information for 
        later methods.
        """
        self.header = fits.getheader(self.filename, ext=0)
        self.science_header = fits.getheader(self.filename, ext=1)
        
        self.wv_units = u.micron
        self.flux_units = u.megajansky / u.steradian
        
        self.reference_wv = self.science_header["CRVAL3"] * self.wv_units 
        self.delta_wv = self.science_header["CDELT3"] * self.wv_units
        
    
    
    def _read_data(self):
        """ 
        Runs a routine which reads relevant data information for 
        later methods.
        """
        self.science_data = fits.getdata(self.filename, ext=1)
        self.error_data = fits.getdata(self.filename, ext=2)
        self.dataquality_flag = fits.getdata(self.filename, ext=3)
    
    
    def high_level_header(self):
        """ 
        Displays a high level overview of the different extension
        headers in the datacube file.
        """
        with fits.open(self.filename) as hdul:
            hdul.info()
    
    
    def index_to_wavelength(self, index):
        """ 
        Returns the wavelength corresponding to a given index.
        
        Arguments
        ---------
        index : int
            The index corresponding to the desired flux measurement.
        
        Returns
        -------
        Quantity : The corresponding wavelength with appropriate units.
        """
        
        if isinstance(index, int):
            if index >= 0 and index < (self.science_data.shape[0]):
                return self.reference_wv + (index + 1) * self.delta_wv
            else:
                raise ValueError(f"Provided index is out of bounds! Recieved {index} while data ranges from {0} to {self.science_data.shape[0] - 1}.")
        else:
            raise TypeError(f"Provided index must be an int! Recieved {type(index)} instead.")
    
    def wavelength_to_index(self, wavelength):
        """ 
        Returns the closest index corresponding to a given wavelength 
        value.
        
        Arguments
        ---------
        wavelength : Quantity
            The wavelength provided in appropriate Astropy units.
        
        Returns
        -------
        int : An integer corresponding to the closest index.
        """
        from astropy.units import Quantity
        
        if isinstance(wavelength, Quantity):
            if wavelength >= self.reference_wv and wavelength < self.reference_wv + self.science_data.shape[0] * self.delta_wv:
                return int(np.floor(((wavelength - self.reference_wv) / self.delta_wv).decompose() - 1))
            else:
                raise ValueError(f"Provided wavelength is out of bounds! Recieved {wavelength} while data ranges from {self.reference_wv} to {self.reference_wv + self.science_data.shape[0] * self.delta_wv}.")
        else:
            raise TypeError(f"Provided wavelength must be an Quantity! Recieved {type(wavelength)} instead.")
