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
#from hydrogenase_processing.prospecpy import ProSpecPy


#Function to Cut the Range into the desired amount
def cut_range_subtraction(raw_spectra:OpusData, raw_wv:OpusData, range_start:int = 3997, range_end:int = 499, SG_poly:int = 3,SG_points:int = 21) -> Tuple[List[AtmFitParams], np.ndarray]:
    """
    Cuts the specified range from raw spectra data and performs atmospheric subtraction.
    
    Args:
        raw_spectra (OpusData): Raw spectral data (FTIR measurements).
        raw_wv (OpusData): Raw water vapor data (water vapor absorption spectra).
        range_start (int, optional): Start wavenumber for the desired range. Defaults to 3997.
        range_end (int, optional): End wavenumber for the desired range. Defaults to 499.
        SG_poly (int, optional): Polynomial order for Savitzky-Golay smoothing. Defaults to 3.
        SG_points (int, optional): Number of points for Savitzky-Golay smoothing. Defaults to 21.
    Returns:
        Tuple[List[AtmFitParams], np.ndarray]: Atmospheric subtraction results.
            -List of atmospheric fitting parameters (e.g. coefficients for water vapor correction.)
            -NumPy array representing the results of atmospheric subtraction.
    Notes: 
        - The function extracts wavenumbers from the raw data and water vapor.
        - It identifies the specified wavenumber range within the raw spectra.
        - Savitzky-Golay smoothing is applied to enhance data quality.
        - The atmospheric subtraction results are returned.
    """
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

    #print(range_start, range_end)
    #print(ind_range_start, ind_range_end)
 

    raw_wavenb_cut = raw_spectra_x[ind_range_start:ind_range_end + 1]
    raw_ab_cut = raw_spectra_abs[ind_range_start:ind_range_end +  1]
    wv_wavenb_cut = wv_x[ind_range_start:ind_range_end + 1]
    wv_ab_cut = raw_wv_abs[ind_range_start:ind_range_end + 1]


    cut_raw_sub_cut_wv = atm_subtraction(raw_wavenb_cut, raw_ab_cut, wv_ab_cut, SG_poly, SG_points)

    return cut_raw_sub_cut_wv

def cut_range_subtraction_multiple_wv(raw_spectra:OpusData, raw_wv:dict, range_start:int = 3997, range_end:int = 499, SG_poly:int = 3,SG_points:int = 21) -> Tuple[List[AtmFitParams], np.ndarray]:
    """
    Cuts the specified range from raw spectra data and performs atmospheric subtraction using multiple water vapor spectra.
    
    Args:
        raw_spectra (OpusData): Raw spectral data (FTIR measurements).
        raw_wv (dict): Dict of all water vapor spectra, with appropriately named keys.
        range_start (int, optional): Start wavenumber for the desired range. Defaults to 3997.
        range_end (int, optional): End wavenumber for the desired range. Defaults to 499.
        SG_poly (int, optional): Polynomial order for Savitzky-Golay smoothing. Defaults to 3.
        SG_points (int, optional): Number of points for Savitzky-Golay smoothing. Defaults to 21.
    Returns:
        Tuple[List[AtmFitParams], np.ndarray]: Atmospheric subtraction results.
            -List of atmospheric fitting parameters (e.g. coefficients for water vapor correction.)
            -NumPy array representing the results of atmospheric subtraction.
    Notes: 
        - The function extracts wavenumbers from the raw data and water vapor, cuttting them into range, and also combining the water vapor columnwise.
        - It identifies the specified wavenumber range within the raw spectra.
        - Savitzky-Golay smoothing is applied to enhance data quality.
        - The atmospheric subtraction results are returned.
    """
    #Extracting the wavenumbers from the raw data and water vapor
    raw_spectra_x = raw_spectra.get_range("AB")
    # the "AB" data can contain more null values at the end (at least 1)
    # so the getting useful data requires slicing the array:
    raw_spectra_abs = raw_spectra["AB"][0:len(raw_spectra_x)]

    #Empty dict of all the absorbances from each vapor spectra
    wv_abs = dict()
    
    #Filling the dict with each absorbance
    for i in raw_wv:
        wv_data_i = raw_wv[i]
        wv_abs[i] = wv_data_i["AB"][0:len(wv_data_i.get_range("AB"))]

    

    #Rounding the wavenumbers to the nearest whole number in order to identify the range of interest.
    whole_num_raw_spectra_x = np.round(raw_spectra_x, 0)

    #Extracting the indices of the range of interest within the NumPy array of the raw spectrum wavenumbers
    ind_range_start = np.where(np.logical_and(whole_num_raw_spectra_x >= range_start - 2, whole_num_raw_spectra_x < range_start + 1))[0][0]

    ind_range_end = np.where(np.logical_and(whole_num_raw_spectra_x >= range_end - 2, whole_num_raw_spectra_x < range_end + 1))[0][0]

    #Empty dict of all cut abs data
    wv_abs_cut = dict() 

    #Filling dict with cut abs data
    for i in wv_abs:
        wv_abs_cut[f'{i}_cut'] = wv_abs[i][ind_range_start:ind_range_end + 1]

    #Combining all cut wv abs data into one columnwise array
    wv_combined_cut_abs = np.column_stack(list(wv_abs_cut.values()))

    #print(range_start, range_end)
    print(ind_range_start, ind_range_end)
 

    raw_wavenb_cut = raw_spectra_x[ind_range_start:ind_range_end + 1]
    raw_ab_cut = raw_spectra_abs[ind_range_start:ind_range_end +  1]


    cut_raw_sub_cut_wv = atm_subtraction(raw_wavenb_cut, raw_ab_cut, wv_combined_cut_abs, SG_poly, SG_points)

    return cut_raw_sub_cut_wv


def cut_range_subtract_multiple_wv(raw_spectra: List, raw_wv:dict, range_start:int = 3997, range_end:int = 499, SG_poly:int = 3, SG_points:int = 21) -> dict[OpusData]:
    """
    Cuts the specified range from raw spectra data and performs atmospheric subtraction over an entire dict of raw spectra.
    
    Args:
        raw_spectra dict[OpusData]:Multiple raw spectra (FTIR measurements).
        raw_wv (OpusData): Multiple water vapor data (water vapor absorption spectra).
        range_start (int, optional): Start wavenumber for the desired range. Defaults to 3997.
        range_end (int, optional): End wavenumber for the desired range. Defaults to 499.
        SG_poly (int, optional): Polynomial order for Savitzky-Golay smoothing. Defaults to 3.
        SG_points (int, optional): Number of points for Savitzky-Golay smoothing. Defaults to 21.
    Returns:
        dict[OpusData]: OpusData objects cut to range with water vapor subtracted

    Notes: 
        - The function extracts wavenumbers from the raw data and water vapor.
        - It identifies the specified wavenumber range within the raw spectra.
        - Savitzky-Golay smoothing is applied to enhance data quality.
        - The atmospheric subtraction results are returned.
    """
    cut_sub_wv_data = dict()

    for prospecpy_obj in raw_spectra:
        prospecpy_obj.cut_range_subtract(raw_wv, range_start, range_end)
        """cut_sub_data = cut_range_subtraction_multiple_wv(prospecpy_obj.get_raw_data(), raw_wv, range_start, range_end)
        #print(cut_sub_data[1])
        prospecpy_obj.cut_atmfitparams_obj = cut_sub_data[0]
        prospecpy_obj.cut_subtracted_data = cut_sub_data[1]
        #cut_sub_wv_data[index] = cut_range_subtraction_multiple_wv(prospecpy_obj.get_raw_data(), raw_wv, range_start, range_end)
        """
    return cut_sub_wv_data