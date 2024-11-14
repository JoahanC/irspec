from cubespec import CubeSpec
from astropy.io import fits
import astropy.units as u


#spec_obj = CubeSpec("./../", "param_files", "IR23128-S_11_grid_param.txt", "input_data/IR23128-S/", redshift=0.044601, fit_dirname="IR23128S_AGN11", mode="AGN")
#spec_obj.perform_grid_extraction()

grid_output_fitsfile = "./../creta_output/extractions/IR23128-S/7/IR23128-S_GridExt_13x13_s1.0as.fits"

fitsob = fits.open(grid_output_fitsfile)
fitsob_header = fits.getheader(grid_output_fitsfile, ext=1)
fitsob_data = fits.getdata(grid_output_fitsfile)
#print(len(fitsob_header.cards))


def construct_header_keys(header):
    cards = header.cards
    
    grid_information = {}
    ra_dec_to_pix = {}
    for card in cards:
        if card[0].isnumeric():
            card_values = card[1].split()
            key_tuple = (card_values[2], card_values[6])
            if key_tuple not in ra_dec_to_pix:
                ra_dec_to_pix[key_tuple] = [card_values[33][:-1], card_values[35][:-1]]
            if key_tuple not in grid_information:
                grid_information[(card_values[33][:-1], card_values[35][:-1])] = [card_values[2] * u.degree, card_values[6] * u.degree, card_values[10], card_values[14]]
    return grid_information, ra_dec_to_pix
            

grid_information, ra_dec_to_pix = construct_header_keys(fitsob_header)
print(grid_information)
#print(repr(fitsob_header))
#print(len(fitsob_data))
#print(fitsob_data[0])