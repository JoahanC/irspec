from cubespec import CubeSpec
#import astropy.units as u

core = "S"
mode = "AGN"
index = 3

#spec_obj = CubeSpec("./../", "param_files", "miri_ir15250_param.txt", "input_data/IR15250/", redshift=0.055206, dirname="IR15250_AGN", asec="03")
#spec_obj = CubeSpec("./../", "param_files", "miri_ir23128n_param.txt", "input_data/IR23128-N/", redshift=0.044601, dirname="IR23128-N_AGN", asec="03")
#spec_obj = CubeSpec("./../", "param_files", "miri_ir23128s_param.txt", "input_data/IR23128-S/", redshift=0.044601, dirname="IR23128-S_AGN", asec="07")

#spec_obj.recall_line()
spec_obj = CubeSpec("./../", "param_files", f"IR23128-{core}_{index}_single_param.txt", f"input_data/IR23128-{core}/", redshift=0.044601, fit_dirname=f"IR23128{core}_{mode}{index}", mode=mode)

#spec_obj.recall_fit()
#spec_obj.perform_fit()
#spec_obj.perform_single_extraction_custom()
#spec_obj.recall_fit()

#spec_obj = CubeSpec("./../", "param_files", "IR23128-S_0_single_param.txt", "input_data/IR23128-S/", redshift=0.044601, fit_dirname="IR23128S_AGN0", mode="AGN")
#spec_obj.perform_single_extraction_custom()
#spec_obj.perform_fit()
spec_obj.recall_line()
