"""Readers for the packaged line-parameter tables in spaxel_fit_params/."""
from importlib.resources import files


def _param_file(name):
    return (files("irspec") / "spaxel_fit_params" / name).open("r")


def read_line_params():
    with _param_file("fitparams.dat") as filey:
        param_data = filey.readlines()
    line_params = {}
    for line in param_data:
        datums = line.split()
        line_params[datums[0]] = [
            int(datums[1]), datums[2], [float(datums[3]), float(datums[4])], 
            [float(datums[5]), float(datums[6])], float(datums[7]), 
            float(datums[8]), float(datums[9]), float(datums[10]), 
            float(datums[11]), float(datums[12]), float(datums[13])]
        if len(datums) == 15:
            line_params[datums[0]].append(float(datums[14]))
    return line_params

def read_line_params2():
    with _param_file("fitparams_n.dat") as filey:
        param_data = filey.readlines()
    line_params = {}
    for line in param_data:
        datums = line.split()
        line_params[datums[0]] = [
            int(datums[1]), datums[2], [float(datums[3]), float(datums[4])], 
            [float(datums[5]), float(datums[6])], float(datums[7]), 
            float(datums[8]), float(datums[9]), float(datums[10]), 
            float(datums[11]), float(datums[12]), float(datums[13])]
        if len(datums) == 15:
            line_params[datums[0]].append(float(datums[14]))
    return line_params

def read_line_params3():
    with _param_file("fitparams_l.dat") as filey:
        param_data = filey.readlines()
    line_params = {}
    for line in param_data:
        datums = line.split()
        line_params[datums[0]] = [
            int(datums[1]), datums[2], [float(datums[3]), float(datums[4])], 
            [float(datums[5]), float(datums[6])], float(datums[7]), 
            float(datums[8]), float(datums[9]), float(datums[10]), 
            float(datums[11]), float(datums[12]), float(datums[13])]
        if len(datums) == 15:
            line_params[datums[0]].append(float(datums[14]))
    return line_params