#Setting Up and Importing the Necessary Packages/Libraries
##Package for reading in Bruker OPUS type files
from brukeropusreader import read_file, OpusData
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
##Package for subtraction algorithm
from scipy.signal import savgol_filter as savitzky_golay
from scipy.optimize import least_squares as leastsq
#Local Files
from .vaporfit import atm_subtraction
from .vaporfit import AtmFitParams
from typing import Tuple, List, SupportsInt, AnyStr


#Function to Cut the Range into the desired amount
#Katelyn initial docstrings
def cut_range_subtraction(raw_spectra:OpusData, raw_wv:OpusData, range_start:int = 3997, range_end:int = 499, SG_poly:int = 3,SG_points:int = 21) -> Tuple[List[AtmFitParams], np.ndarray]:
    #Extracting the wavenumbers from the raw data and water vapor
    raw_spectra_x = raw_spectra.get_range("AB")
    # the "AB" data can contain more null values at the end (at least 1)
    # so the getting useful data requires slicing the array:
    raw_spectra_abs = raw_spectra["AB"][0:len(raw_spectra_x)]

    wv_x = raw_wv.get_range("AB")
    raw_wv_abs = raw_wv["AB"][0:len(wv_x)]

    #Rounding the wavenumbers to the nearest whole number in order to identify the range of interest.
    whole_num_raw_spectra_x = np.round(raw_spectra_x, 0)

    #Extracting the indices of the range of interest within the NumPy array of the raw spectrum wavenumbers
    ind_range_start = np.where(np.logical_and(whole_num_raw_spectra_x >= range_start - 2, whole_num_raw_spectra_x < range_start + 1))[0][0]
    ind_range_end = np.where(np.logical_and(whole_num_raw_spectra_x >= range_end - 2, whole_num_raw_spectra_x < range_end + 1))[0][0]

    print(range_start, range_end)
    print(ind_range_start, ind_range_end)
 

    raw_wavenb_cut = raw_spectra_x[ind_range_start:ind_range_end + 1]
    raw_ab_cut = raw_spectra_abs[ind_range_start:ind_range_end +  1]
    wv_wavenb_cut = wv_x[ind_range_start:ind_range_end + 1]
    wv_ab_cut = raw_wv_abs[ind_range_start:ind_range_end + 1]


    cut_raw_sub_cut_wv = atm_subtraction(raw_wavenb_cut, raw_ab_cut, wv_ab_cut, SG_poly, SG_points)

    return cut_raw_sub_cut_wv


