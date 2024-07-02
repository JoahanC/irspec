from datacube import Datacube
from cubespec import CubeSpec
import astropy.units as u

#0.044601
#0.055206
spec_obj = CubeSpec("./../", "param_files", "miri_ir23128s_param.txt", "input_data/IR23128-S/", redshift=0.044601)
#spec_obj.load_line_params()
#spec_obj.recall_fit()
spec_obj.line_cutouts()
