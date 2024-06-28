from datacube import Datacube
from cubespec import CubeSpec
import astropy.units as u

spec_obj = CubeSpec("./../", "param_files", "miri_ir23128s_param.txt", "input_data/IR23128-S/")
spec_obj.rewrite_spec_csv()

