import glob 
import os
import matplotlib.pyplot as plt
from regions import PolygonSkyRegion, PixCoord
import astropy.units as u
from astropy.coordinates import SkyCoord
from irspec.datacube import Datacube
from irspec.cubespec import CubeSpec


test_path = glob.glob(f"./irspec/data/ch4-long_s3d.fits")[0]
testcube = Datacube(test_path, redshift=0.044601, verbose=True)
spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/data/", "./irspec/param_files/", 
                    "ir23128s_single1S_param.txt", "./../../creta_extractions/single1S",
                    "./../../cafe_output/single1S_AGN", 0.044601, testcube, mode="AGN")

#spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/data/", "./irspec/param_files/", 
#                    "ir23128s_single1S_param.txt", "./../../creta_extractions/single1S",
#                    "./../../cafe_output/single1S_SB_and_AGN", 0.044601, testcube, mode="SB_and_AGN")
spec_obj.recall_line()
#spec_obj.perform_single_extraction()
#spec_obj.perform_fit()
#spec_obj.recall_fit()
#spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/data/", "./irspec/param_files/", 
#                    "ir23128s_single1S_param.txt", "./../../creta_extractions/single1S",
#                    "./../../cafe_output/single1S_SB", 0.044601, testcube, mode="SB")
#spec_obj.perform_single_extraction()
#spec_obj.perform_fit()
def create_directory_if_not_exists(path):
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Directory '{path}' created.")
    else:
        print(f"Directory '{path}' already exists.")

"""for i in range(50):
    
    create_directory_if_not_exists(f"./../../creta_extractions/manual/{i}_single")
    create_directory_if_not_exists(f"./../../cafe_output/manual/{i}_single")
    spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/data/", "./irspec/param_files/manual_grid/", 
                        f"{i}_single_param.txt", f"./../../creta_extractions/manual/{i}_single",
                        f"./../../cafe_output/manual/{i}_single", 0.044601, testcube, mode="SB_and_AGN")
    spec_obj.perform_single_extraction()
    try:
        spec_obj.perform_fit()
    except:
        print("nan_ap")"""

#spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/data/", "./irspec/param_files/", 
#                    "ir23128s_grid3_param.txt", "./../../creta_extractions/fe2ar2",
#                    "./../../cafe_output/single1S_AGN", 0.044601, testcube, mode="AGN")
#spec_obj.perform_grid_extraction()
