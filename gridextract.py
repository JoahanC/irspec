from cubespec import CubeSpec
from astropy.io import fits
#import astropy.units as u


#spec_obj = CubeSpec("./../", "param_files", "IR23128-S_11_grid_param.txt", "input_data/IR23128-S/", redshift=0.044601, fit_dirname="IR23128S_AGN11", mode="AGN")
#spec_obj.perform_grid_extraction()

grid_output_fitsfile = "./../creta_output/extractions/IR23128-S/7/IR23128-S_GridExt_13x13_s1.0as.fits"

fitsob = fits.open(grid_output_fitsfile)
fitsob_header = fits.getheader(grid_output_fitsfile, ext=1)
fitsob_data = fits.getdata(grid_output_fitsfile)
#print(len(fitsob_header.cards))

def construct_header_keys(header):
    cards = header.cards
    
    for card in cards:
        if card[0].isnumeric():
            print(card[1])

construct_header_keys(fitsob_header)
#print(repr(fitsob_header))
#print(len(fitsob_data))
#print(fitsob_data[0])