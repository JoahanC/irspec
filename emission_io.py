def read_line_params():
    with open("./fitparams.dat", 'r') as filey:
        param_data = filey.readlines()
    line_params = {}
    for line in param_data:
        datums = line.split()
        line_params[datums[0]] = [int(datums[1]), datums[2], [float(datums[3]), float(datums[4])], [float(datums[5]), float(datums[6])], float(datums[7]), [float(datums[8]), float(datums[9]), float(datums[10]), float(datums[11])]]
    return line_params
