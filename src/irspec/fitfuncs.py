""" 
This file contains all of the fitting functions used in CubeSpec
"""
import numpy as np


def OneDPolynomial(x, a, b):
        return a * x + b 
    
def TwoDPolynomial(x, a, b, c):
    return a * x ** 2 + b * x + c 

def ThreeDPolynomial(x, a, b, c, d):
    return a * x ** 3 + b * x ** 2 + c * x + d

def OneDGaussian(x, A, x0, sigma):
    return A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

def PowerLaw(x, A, k, c):
    return A * x ** k + c