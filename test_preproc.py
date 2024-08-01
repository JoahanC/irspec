""" 
Testing script for diagnosing the performance of the modified 
cube_preproc class.
"""
import numpy as np
from astropy.coordinates import SkyCoord

from cube_preproc import cube_preproc
cubepath = "./../input_data/IR23128-N/Level3_ch1-short_s3d.fits"
testproc = cube_preproc()
res = testproc.getFITSData(cubepath)
ra = "23h15m46.7878s"
dec = "-59d03m09.770s"
pos1 = SkyCoord(ra, dec, frame='icrs')
testproc.plotGrid(cubepath, pos1.ra.deg, pos1.dec.deg, 0.5, 10, 10, -1, "test", "testy")
