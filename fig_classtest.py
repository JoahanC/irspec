from localfit import LocalFit
from cubespec import CubeSpec

spec_obj = CubeSpec("./../", "param_files", "IR23128-S_0_single_param.txt", "input_data/IR23128-S/", redshift=0.044601, fit_dirname="IR23128S_AGN0", mode="AGN")
data = spec_obj.recall_data()


"""figob = LocalFit(data["wave"], data["flux"], [15, 16], 15.553, "[NeIII]")
figob.main_fit(1, 3)
figob.render_fit()"""


figob = LocalFit(data["wave"], data["flux"], [18.4, 19], 18.713, "[SIII]")
figob.main_fit(1, 1)
figob.render_fit()

"""figob = LocalFit(data["wave"], data["flux"], [12.5, 13.1], 12.813, "[NeII]")
figob.main_fit(2, 2)
figob.render_fit()"""

"""figob = LocalFit(data["wave"], data["flux"], [10.27, 10.84], 10.507, "[SIV]")
figob.main_fit(2, 2)
figob.render_fit()"""