""" 
This is a testing file to be used for diagnosing the performance of 
the Datacube class.
"""
from datacube import Datacube

filepath = "./../input_data/IR15250/Level3_ch1-long_s3d.fits"

testcube = Datacube(filepath)
