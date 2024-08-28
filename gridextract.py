from cubespec import CubeSpec
#import astropy.units as u


spec_obj = CubeSpec("./../", "param_files", "IR23128-S_11_grid_param.txt", "input_data/IR23128-S/", redshift=0.044601, fit_dirname="IR23128S_AGN11", mode="AGN")
spec_obj.perform_grid_extraction()