import glob
import warnings
from irspec.datacube import Datacube
from irspec.spaxel_fit import Spaxelcube
from irspec.emission_io import read_line_params, read_line_params2
from irspec.sort_dynamic import sortdata
from astropy.io import fits
from irspec.cubespec import CubeSpec
from astropy.utils.exceptions import AstropyDeprecationWarning


warnings.filterwarnings('ignore', category=AstropyDeprecationWarning, append=True)
# What instrument to be used
#instrument = "NIRSpec"
instrument = "MRS"

# Read in line parameters
line_dict = read_line_params()

# What lines to fit
#lines = ["Pfund65", "[MgVII]", "[MgV]", "[NeVI]", "[NeV]", "[OIV]", "[SIII]", "[SIV]", "[NeII]", "[NeIII]", "[ArIII]", "[ArII]", "[NeV]_14", "[FeII]", "[H_2_S_1]", "[H_2_S_2]", "[H_2_S_3]", "[H_2_S_5]"]
#lines = ["[MgV]", "[NeVI]", "[NeV]", "[OIV]", "[SIII]", "[SIV]", "[NeII]", "[NeIII]", "[ArIII]", "[ArII]", "[NeV]_14", "[FeII]", "[H_2_S_1]", "[H_2_S_2]", "[H_2_S_3]", "[H_2_S_5]"]
lines = ["[H_2_S_1]", "[H_2_S_2]"]
lines = ["NeIII"]
#lines = ["[KVII]", "[CaIV]", "[AlVI]", "Pfund85", "Brackett54", "[MgIV]", "[ArVI]"]
#lines = ["Brackett54"]
#lines = ["[NeIII]"]#["[NeII]", "[NeIII]"]

# North or south data paths
data_source = ["south"]
source = "south"
diagnose=False

if instrument == "MRS":
    for source in data_source:
        for line in lines:
            channel = line_dict[line][0]
            subchannel = line_dict[line][1]
            print(channel)
            print(subchannel)
            if source == "north":
                print(source)
                print(line)
                test_path = glob.glob(f"./irspec/north_data/ch{channel}-{subchannel}_s3d.fits")[0]
                testcube = Datacube(test_path, redshift=0.044601, verbose=True)
                #print(testcube.area_sr, testcube._flux_unit)
                #test_path = "./../../creta_extractions/fe2ar2/IR23128-S_GridExt_50x50_s0.17000000178813898as_cube.fits"
                #testcube = Datacube(test_path, redshift=0.044601, verbose=True)
                #print(testcube.area_sr, testcube._flux_unit)

                #print(repr(testcube.science_header))
                #testcube.display_dq(300)
                testfit = Spaxelcube(testcube, line, f"./../../spaxs/north/{line}/")
                file_stem = f"./../../spaxs/north/{line}/"
                #testfit.load_fit(f"./../../spaxs/north/{line}/twogaussian_raw.dat")
                testfit.two_gaussian_fit()
                testfit.render_multicomponent_plot(file_stem, savefig=True)
                testfit.render_totflux_plot(file_stem, savefig=True)
                testfit.render_rel_vel_plot("G1CEN", file_stem, savefig=True)
                testfit.render_rel_vel_plot("G2CEN", file_stem, savefig=True)
                testfit.render_sigma_plot("G1SIGMA", file_stem, savefig=True)
                testfit.render_sigma_plot("G2SIGMA", file_stem, savefig=True)
                #testfit.render_spaxel_fit(27, 20)
            if source == "south":
                print(source)
                print(line)
                test_path = glob.glob(f"./irspec/data/ch{channel}-{subchannel}_s3d.fits")[0]
                testcube = Datacube(test_path, redshift=0.044601, verbose=True)
                #print(testcube.area_sr, testcube._flux_unit)
                #test_path = "./../../creta_extractions/fe2ar2/IR23128-S_GridExt_50x50_s0.17000000178813898as_cube.fits"
                #testcube = Datacube(test_path, redshift=0.044601, verbose=True)
                #print(testcube.area_sr, testcube._flux_unit)

                #print(repr(testcube.science_header))
                #testcube.display_dq(300)
                if diagnose:
                    testfit = Spaxelcube(testcube, line, f"./../../spaxs/{line}/", mode="spaxel", test_spaxel=(33,18))
                    #testfit.two_gaussian_fit(diagnose=True)
                else:
                    testfit = Spaxelcube(testcube, line, f"./../../spaxs/{line}/")
                    #testfit.two_gaussian_fit(diagnose=False)
                    testfit.load_fit(f"./../../spaxs/{line}/twogaussian_raw.dat")
                    file_stem = f"./../../spaxs/{line}/"
                    testfit.render_multicomponent_plot(file_stem, savefig=True)
                    testfit.render_totflux_plot(file_stem, savefig=True)
                    testfit.render_rel_vel_plot("G1CEN", file_stem, savefig=True)
                    testfit.render_rel_vel_plot("G2CEN", file_stem, savefig=True)
                    testfit.render_sigma_plot("G1SIGMA", file_stem, savefig=True)
                    testfit.render_sigma_plot("G2SIGMA", file_stem, savefig=True)

#lines = ["Brackett54", "[MgIV]", "[AlVI]"]
#lines = ["[MgIV]", "Brackett54", "[AlVI]", "Pfund85"]
#lines = ["Brackett54"]
if instrument == "NIRSpec":
    for source in data_source:
        for line in lines:
            channel = line_dict[line][0]
            subchannel = line_dict[line][1]
            print(channel)
            print(subchannel)
            #line = "[FeII]"
            if source == "north":
                print(source)
                print(line)
                test_path = glob.glob(f"./irspec/north_data/ch{channel}-{subchannel}_s3d.fits")[0]
                testcube = Datacube(test_path, redshift=0.044601, verbose=True)
                #print(testcube.area_sr, testcube._flux_unit)
                #test_path = "./../../creta_extractions/fe2ar2/IR23128-S_GridExt_50x50_s0.17000000178813898as_cube.fits"
                #testcube = Datacube(test_path, redshift=0.044601, verbose=True)
                #print(testcube.area_sr, testcube._flux_unit)

                #print(repr(testcube.science_header))
                #testcube.display_dq(300)
                testfit = Spaxelcube(testcube, line, f"./../../spaxs/north/{line}/")
                file_stem = f"./../../spaxs/north/{line}/"
                testfit.load_fit(f"./../../spaxs/north/{line}/twogaussian_raw.dat")
                #testfit.two_gaussian_fit()
                testfit.render_multicomponent_plot(file_stem, savefig=True)
                testfit.render_totflux_plot(file_stem, savefig=True)
                testfit.render_rel_vel_plot("G1CEN", file_stem, savefig=True)
                testfit.render_rel_vel_plot("G2CEN", file_stem, savefig=True)
                testfit.render_sigma_plot("G1SIGMA", file_stem, savefig=True)
                testfit.render_sigma_plot("G2SIGMA", file_stem, savefig=True)
                #testfit.render_spaxel_fit(27, 20)
            if source == "south":
                print(source)
                print(line)
                test_path = glob.glob("./irspec/data/g395h-f290lp_s3d.fits")[0]
                testcube = Datacube(test_path, redshift=0.044601, verbose=True)
                #print(testcube.area_sr, testcube._flux_unit)
                #test_path = "./../../creta_extractions/fe2ar2/IR23128-S_GridExt_50x50_s0.17000000178813898as_cube.fits"
                #testcube = Datacube(test_path, redshift=0.044601, verbose=True)
                #print(testcube.area_sr, testcube._flux_unit)

                #print(repr(testcube.science_header))
                #testcube.display_dq(300)
                testfit = Spaxelcube(testcube, line, f"./../../spaxs/{line}/")
                file_stem = f"./../../spaxs/{line}/"
                
                #testfit.two_gaussian_fit()
                testfit.load_fit(f"./../../spaxs/{line}/twogaussian_raw.dat")
                sortdata("CEN", line, line_dict[line][4])
                testfit.load_fit(f"./../../spaxs/{line}/CEN_sort.dat")
                
                
                testfit.render_multicomponent_plot(file_stem, savefig=True)
                testfit.render_totflux_plot(file_stem, savefig=True)
                testfit.render_rel_vel_plot("G1CEN", file_stem, savefig=True)
                testfit.render_rel_vel_plot("G2CEN", file_stem, savefig=True)
                testfit.render_sigma_plot("G1SIGMA", file_stem, savefig=True)
                testfit.render_sigma_plot("G2SIGMA", file_stem, savefig=True)
                #testfit.render_spaxel_fit(50, 38)
                #print(testcube.wv_to_vel(3.747, 3.741))

    """for line in lines:
        #line = "[Brackett54]"
        test_path = glob.glob("./irspec/data/g395h-f290lp_s3d.fits")[0]
        testcube = Datacube(test_path, redshift=0.044601, verbose=True)
        testfit = Spaxelcube(testcube, line, f"./../../spaxs/{line}/")
        file_stem = f"./../../spaxs/{line}/"        
        testfit.load_fit(f"./../../spaxs/{line}/twogaussian_raw.dat")
        #testfit.two_gaussian_fit()
        #testfit.render_multicomponent_plot(file_stem, savefig=True)
        #testfit.render_totflux_plot(file_stem, savefig=False)
        testfit.render_rel_vel_plot("G1CEN", file_stem, savefig=False)
        #testfit.render_rel_vel_plot("G2CEN", file_stem, savefig=True)
        #testfit.render_sigma_plot("G1SIGMA", file_stem, savefig=True)
        #testfit.render_sigma_plot("G2SIGMA", file_stem, savefig=True)
        testfit.render_spaxel_fit(39, 47)
        print(testcube.wv_to_vel(3.749, 3.7406))#"""

test_path = glob.glob(f"./irspec/data/ch4-long_s3d.fits")[0]
testcube = Datacube(test_path, redshift=0.044601, verbose=True)
spec_obj = CubeSpec("/Users/jcj/Documents/research/goals/irspec/src/irspec/data/", "./irspec/param_files/", 
                    "ir23128s_grid2_param.txt", "./../../creta_extractions/ne3ne2/",
                    "./../../cafe_output/grid1_AGN/", 0.044601, testcube, mode="AGN")
#spec_obj.perform_grid_extraction()

"""lines = ["[FeII]", "[ArII]"]
for line in lines:
    channel = line_dict[line][0]
    subchannel = line_dict[line][1]
    test_path = "./../../creta_extractions/fe2ar2/IR23128-S_GridExt_20x20_s0.25as_cube.fits"

    with fits.open(test_path, mode='update') as filehandle:
        filehandle[1].header['CTYPE2'] = "DEC--TAN"
    #print(repr(hdu.header))
    testcube = Datacube(test_path, redshift=0.044601, verbose=True)
    print(testcube.wvs)
    #print(testcube.area_sr, testcube._flux_unit)

    #print(repr(testcube.science_header))
    #testcube.display_dq(300)
    testfit = Spaxelcube(testcube, line, f"./../../spaxs/test/{line}/")
    file_stem = f"./../../spaxs/test/{line}/"
    testfit.load_fit(f"./../../spaxs/test/{line}/twogaussian_raw.dat")
    #testfit.two_gaussian_fit()
    testfit.render_multicomponent_plot(file_stem, savefig=True)
    testfit.render_totflux_plot(file_stem, savefig=True)
    testfit.render_rel_vel_plot("G1CEN", file_stem, savefig=True)
    testfit.render_rel_vel_plot("G2CEN", file_stem, savefig=True)
    testfit.render_sigma_plot("G1SIGMA", file_stem, savefig=True)
    testfit.render_sigma_plot("G2SIGMA", file_stem, savefig=True)
    #testfit.render_spaxel_fit(3, 3)
    #print(testcube.wvs[-1])"""

#lines = ["[NeV]_14", "[NeIII]", "[NeII]"]
#lines = ["[ArII]"]#, "[FeII]"]
if instrument == "ratios":
    lines = ["[NeII]"]
    files = []
    for line in lines:
        channel = line_dict[line][0]
        subchannel = line_dict[line][1]
        test_path = "./../../creta_extractions/ne3ne2/IR23128-S_GridExt_40x40_s0.273as_cube.fits"
        #test_path = "./../../creta_extractions/fe2ar2/IR23128-S_GridExt_35x35_s0.196as_cube.fits"
        #test_path = glob.glob(f"./irspec/data/ch{channel}-{subchannel}_s3d.fits")[0]

        with fits.open(test_path, mode='update') as filehandle:
            filehandle[1].header['CTYPE2'] = "DEC--TAN"
        #print(repr(hdu.header))
        testcube = Datacube(test_path, redshift=0.044601, verbose=True)
        print(testcube.wvs)
        #print(testcube.area_sr, testcube._flux_unit)

        #print(repr(testcube.science_header))
        #testcube.display_dq(300)
        testfit = Spaxelcube(testcube, line, f"./../../spaxs/test/{line}/")
        file_stem = f"./../../spaxs/test/{line}/"
        #testfit.load_fit(f"./../../spaxs/test/{line}/twogaussian_raw.dat")
        
        testfit.two_gaussian_fit()
        testfit.render_multicomponent_plot(file_stem, savefig=True)
        testfit.render_totflux_plot(file_stem, savefig=True)
        testfit.render_rel_vel_plot("G1CEN", file_stem, savefig=True)
        testfit.render_rel_vel_plot("G2CEN", file_stem, savefig=True)
        testfit.render_sigma_plot("G1SIGMA", file_stem, savefig=True)
        testfit.render_sigma_plot("G2SIGMA", file_stem, savefig=True)
        #testfit.render_spaxel_fit(30, 15)
        #print(testcube.wvs[-1])