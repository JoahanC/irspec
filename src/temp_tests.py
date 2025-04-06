import glob
from irspec.datacube import Datacube
from irspec.spaxel_fit import SpaxelFit

instrument = "MRS"

if instrument == "MRS":
    channel = "2"
    subchannel = "long"
    test_path = glob.glob(f"./irspec/data/ch{channel}-{subchannel}_s3d.fits")[0]
    testcube = Datacube(test_path, redshift=0.044601, verbose=True)

    #print(repr(testcube.science_header))
    #testcube.display_dq(300)
    #print(testcube.vel_to_sigma(12.83, 1000))
    testfit = SpaxelFit(testcube, "[SIV]", "./../../spaxs/[SIV]/")
    testfit.load_fit("./../../spaxs/[SIV]/twogaussian_raw.dat")
    #testfit.two_gaussian_fit()
    testfit.render_multicomponent_plot()
    testfit.render_totflux_plot()
    testfit.render_rel_vel_plot("G2CEN")
    #testfit.render_spaxel_fit(35, 16)

if instrument == "NIRSpec":
    pass

