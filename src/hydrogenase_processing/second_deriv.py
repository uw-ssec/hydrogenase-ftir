#Setting Up and Importing the Necessary Packages/Libraries
##Package for reading in Bruker OPUS type files
from brukeropusreader import read_file, OpusData
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
##Package for subtraction algorithm
from scipy.signal import savgol_filter as savitzky_golay
from scipy.optimize import least_squares as leastsq
from scipy.interpolate import UnivariateSpline
from scipy.interpolate._fitpack2 import InterpolatedUnivariateSpline
#Local Files
from .vaporfit import atm_subtraction
from .vaporfit import AtmFitParams
from typing import Tuple, List, SupportsInt, AnyStr


#flip element order function: input array and name for the inverted array
def flip_order(x):
    arr_invert = x[::-1]
    return arr_invert

#This function returns a spline representation of the 2nd deriv of the data passed in and a numpy array with that spline expressed over a range, for passing into find peaks.
def second_deriv(cut_sub_data:Tuple[List[AtmFitParams], np.ndarray], show_plots = True):
    ## Extracting wavenumbers and y absorbance, and plotting it
    #wavenumbers
    x_wavenb = cut_sub_data[0][0].wavenb
    #absorbance y values
    y_corr_abs = cut_sub_data[0][0].sub_spectrum
    
    ## Taking Second Derivative, using GRADIENT
    dydx = np.gradient(y_corr_abs, x_wavenb)
    d2ydx2 = np.gradient(dydx, x_wavenb)
    
    #Inverting Data, as UnivariateSpline only takes strictly increasing data
    x_wavenb_invert = flip_order(x_wavenb)
    d2ydx2_invert = flip_order(d2ydx2)
    x_range = np.linspace(x_wavenb[0],x_wavenb[-1],1000)

    d2ydx2_spl = UnivariateSpline(x_wavenb_invert,d2ydx2_invert,s=0,k=3)

    
    
    if show_plots:
        plt.plot(x_wavenb, y_corr_abs)
        plt.title("Plotting Cut and Subtracted Data")
        plt.show()

        plt.plot(x_wavenb, d2ydx2)
        plt.title("Plotting Second Derivative Calculated from Gradient")
        plt.show()

    spline_over_range = d2ydx2_spl(x_range)

    return d2ydx2_spl, spline_over_range, x_range





