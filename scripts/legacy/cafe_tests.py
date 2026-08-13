import glob 
import numpy as np
from irspec.datacube import Datacube
from irspec.cubespec import CubeSpec
import astropy.units as u
import astropy.constants as const

test = False
refit_north = False
refit_south_outflow = False
refit_grid = True
view = False

if test: 
    test_path = glob.glob(f"./irspec/north_data/ch1-short_s3d.fits")[0]
    testcube = Datacube(test_path, redshift=0.044601, verbose=True)
    spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/north_data/", "./irspec/param_files/", 
                        "ir23128s_single2N_param.txt", "./../../creta_extractions/test/",
                        "./../../cafe_output/test/", 0.044601, testcube, mode="SB")
    spec_obj.perform_single_extraction()
    spec_obj.perform_fit()

if refit_north:
    test_path = glob.glob(f"./irspec/north_data/ch1-short_s3d.fits")[0]
    testcube = Datacube(test_path, redshift=0.044601, verbose=True)
    spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/north_data/", "./irspec/param_files/", 
                        "ir23128s_single1N_param.txt", "./../../creta_extractions/single1N/",
                        "./../../cafe_output/single1N_SB/", 0.044601, testcube, mode="SB")
    #spec_obj.perform_single_extraction()
    #spec_obj.perform_fit()
    spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/north_data/", "./irspec/param_files/", 
                        "ir23128s_single2N_param.txt", "./../../creta_extractions/single2N/",
                        "./../../cafe_output/single2N_SB/", 0.044601, testcube, mode="SB")
    spec_obj.perform_single_extraction()
    spec_obj.perform_fit()
if refit_south_outflow:
    test_path = glob.glob(f"./irspec/north_data/ch2-short_s3d.fits")[0]
    testcube = Datacube(test_path, redshift=0.044601, verbose=True)
    spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/data/", "./irspec/param_files/", 
                        "ir23128s_single4S_param.txt", "./../../creta_extractions/single4S/",
                        "./../../cafe_output/single4S_SB/", 0.044601, testcube, mode="SB")
    #spec_obj.perform_single_extraction()
    #spec_obj.perform_fit()
    spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/data/", "./irspec/param_files/", 
                        "ir23128s_single5S_param.txt", "./../../creta_extractions/single5S/",
                        "./../../cafe_output/single5S_SB/", 0.044601, testcube, mode="SB")
    spec_obj.perform_single_extraction()
    spec_obj.perform_fit()
    spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/data/", "./irspec/param_files/", 
                        "ir23128s_single4S_param.txt", "./../../creta_extractions/single4S/",
                        "./../../cafe_output/single4S_AGN/", 0.044601, testcube, mode="AGN")
    #spec_obj.perform_fit()
    spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/data/", "./irspec/param_files/", 
                        "ir23128s_single5S_param.txt", "./../../creta_extractions/single5S/",
                        "./../../cafe_output/single5S_AGN/", 0.044601, testcube, mode="AGN")
    spec_obj.perform_fit()
    spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/data/", "./irspec/param_files/", 
                        "ir23128s_single4S_param.txt", "./../../creta_extractions/single4S/",
                        "./../../cafe_output/single4S_SB_and_AGN/", 0.044601, testcube, mode="SB_and_AGN")
    #spec_obj.perform_fit()
    spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/data/", "./irspec/param_files/", 
                        "ir23128s_single5S_param.txt", "./../../creta_extractions/single5S/",
                        "./../../cafe_output/single5S_SB_and_AGN/", 0.044601, testcube, mode="SB_and_AGN")
    spec_obj.perform_fit()
if refit_grid == True:
    test_path = glob.glob(f"./irspec/data/ch1-short_s3d.fits")[0]
    testcube = Datacube(test_path, redshift=0.044601, verbose=True)
    #spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/data/", "./irspec/param_files/", 
    #                    "ir23128s_grid3_param.txt", "./../../creta_extractions/fe2ar2/",
    #                    "./../../cafe_output/grid1_AGN/", 0.044601, testcube, mode="AGN")
    
    #spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/data/", "./irspec/param_files/", 
    #                    "ir23128s_grid2_param.txt", "./../../creta_extractions/ne3ne2/",
    #                    "./../../cafe_output/grid1_AGN/", 0.044601, testcube, mode="AGN")
    spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/data/", "./irspec/param_files/", 
                        "ir23128s_grid4_param.txt", "./../../creta_extractions/pah/",
                        "./../../cafe_output/grid1_AGN/", 0.044601, testcube, mode="AGN")
    spec_obj.perform_grid_extraction()
if view:
    test_path = glob.glob(f"./irspec/data/ch1-short_s3d.fits")[0]
    testcube = Datacube(test_path, redshift=0.044601, verbose=True)
    spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/data/", "./irspec/param_files/", 
                        "ir23128s_single1S_param.txt", "./../../creta_extractions/single1S/",
                        "./../../cafe_output/single1S_AGN/", 0.044601, testcube, mode="AGN")
    
    '''test_path = glob.glob(f"./irspec/data/ch1-short_s3d.fits")[0]
    testcube = Datacube(test_path, redshift=0.044601, verbose=True)
    spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/data/", "./irspec/param_files/", 
                        "ir23128s_single3S_param.txt", "./../../creta_extractions/single3S/",
                        "./../../cafe_output/single3S_AGN/", 0.044601, testcube, mode="AGN")'''
    
    #asdf = spec_obj.open_asdf()
    #for key in asdf.keys():
    #    print(key, asdf[key].keys())
    #print(asdf["cafefit"]["cont_profs"]['kAbs'])
    #print(asdf["cafefit"]["cont_profs"]['kAbs'])
    #spec_obj.recall_fit()
    #spec_obj.recall_line()
    lines = ["[ArII]", "[ArIII]", "[NeII]", "[NeIII]", "[SIII]", "[SIV]", "[NeV]", "[NeV]_14", "[NeVI]", "[OIV]", "[MgIV]", "[MgV]"]
    amps, cens, sigs = spec_obj.refit_lines(lines)
    
    for idx, amp in enumerate(amps):
        #print(amp.unit)
        #print(cens)
        #gamma = np.abs(sigs[idx]) * 2.355 / cens[idx]
        #flux = (1 / np.sqrt(np.pi * np.log(2)) * (np.pi * const.c.to('micron/s') / 2) * (amps[idx].value * u.Jy * gamma / (cens[idx].value * u.micron))).to(u.watt / u.meter ** 2)
        #flux = (1 / np.sqrt(np.pi * np.log(2)) * (np.pi * const.c.to('micron/s') / 2) * (self.fitparams["G1AMP"][idx] * u.Jy * g1_gamma / (self.fitparams["G1CEN"][idx] * u.micron))).to(u.watt / u.meter ** 2)
        print(f"{lines[idx]} flux: {amp.value}")

#print(testcube.pixel_scale)
#spec_obj.perform_single_extraction()
#spec_obj.perform_grid_extraction()
#path = "./../../creta_extractions/grid1/*.fits"
#grid_file = glob.glob(path)[1]
#print(type(spec_obj.c))
#spec_obj.read_grid_extraction(grid_file)
#spec_obj.perform_fit_cube()
#lines = ["Pfund85", "Brackett54"]
#lines = ["[ArII]", "[ArIII]", "[NeII]", "[NeIII]", "[SIII]", "[SIV]"]
lines = ["[ArII]", "[ArIII]", "[NeII]", "[NeIII]", "[SIII]", "[SIV]", "[NeV]", "[NeV]_14", "[NeVI]", "[OIV]", "[MgIV]", "[MgV]"]
#lines = ["[H_2_S_1]", "[H_2_S_2]", "[H_2_S_3]", "[H_2_S_5]"]
#lines = ["[NeV]", "[NeV]_14", "[NeVI]", "[OIV]", "[MgIV]", "[MgV]"]
#lines = ["Pfund85", "Brackett54","[ArII]", "[ArIII]", "[NeII]", "[NeIII]", "[SIII]", "[SIV]","[H_2_S_1]", "[H_2_S_2]", "[H_2_S_3]", "[H_2_S_5]","[NeV]", "[NeV]_14", "[NeVI]", "[OIV]", "[MgIV]", "[MgV]", "[AlVI]"]
#lines = ["[KVII]", "[CaIV]", "[AlVI]", "[MgIV]"]
#spec_obj.refit_lines(lines)
#spec_obj.recall_line()