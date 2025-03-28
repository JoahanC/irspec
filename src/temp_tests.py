import glob
from irspec.datacube import Datacube
from irspec.spaxel_fit import SpaxelFit


test_path = glob.glob("./irspec/data/*.fits")[0]
testcube = Datacube(test_path, redshift=0.044601, verbose=True)
#print(repr(testcube.science_header))
#testcube.display_dq(300)
testfit = SpaxelFit(testcube, "[NeII]", "./../../spaxs/")
testfit.load_fit("./../../spaxs/twogaussian_raw.dat")
#testfit.two_gaussian_fit()
#testfit.render_multicomponent_plot()
#testfit.render_totflux_plot()
#testfit.render_rel_vel_plot("G2CEN")
testfit.render_spaxel_fit(22, 19)