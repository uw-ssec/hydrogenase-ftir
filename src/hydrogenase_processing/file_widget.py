#Setting Up and Importing the Necessary Packages/Libraries
##Package for reading in Bruker OPUS type files
from brukeropusreader import read_file
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import pathlib
import numpy as np
import pandas as pd
#Local Functions
from hydrogenase_processing.cut_range import cut_range_subtraction, cut_range_subtraction_multiple_wv
from hydrogenase_processing.second_deriv import second_deriv
from hydrogenase_processing.anchor_points import interact
from hydrogenase_processing.anchor_points import get_peaks, get_start_end_anchorpoints, get_all_anchor_points
from hydrogenase_processing.vaporfit import atm_subtraction
from hydrogenase_processing.vaporfit import AtmFitParams
from scipy.signal import find_peaks, peak_widths
from hydrogenase_processing.anchor_points import get_peaks, get_start_end_anchorpoints, get_all_anchor_points, baseline_spline, get_peaks_absorbance, baseline_correction, get_peak_baseline_absorbance, plot_baseline_data
from hydrogenase_processing.peak_fit import gaussian, peak_fit, lorentzian

import ipywidgets as widgets

## pipeline automated up to the adj_factor and threshold
def auto_path(path_to_water_vapor_data, path_to_output_plots, path_to_data, run_code, path_to_excel, sheet_name, sample_raw_name, sample_second_deriv_name):
    
    path_to_water_vapor_data = pathlib.Path(path_to_water_vapor_data)
    path_to_output_plots_= pathlib.Path(path_to_output_plots)
    path_to_data = pathlib.Path(path_to_data)

    #Pulling in all pD6 sample data
    raw_files = list(path_to_data.iterdir())
    raw_files.sort()

    #Initializing dict of raw spectra files from the file system
    raw_data = dict()

    #Populating the raw_test_data dict with all the read in raw opus files
    ##Using the last 5 characters, as they are the uniquely identifying portions of each of the file names
        
    for i in raw_files:
        if not i.name.startswith('.DS_Store'):
            raw_data[f'{run_code}_{i.name[-5:len(i.name)]}'] = read_file(i)

    #print(raw_data.keys())

    #Pulling in all wv data
    water_vapor_files = list(path_to_water_vapor_data.iterdir())
    water_vapor_files.sort()
    #Initializing dict of wv_files from the file system
    water_vapor_data = dict()

    #Populating the water_vapor_data dict with all the read in wv opus files
    #making sure names(keys) are distinct by subscripting
    for i in (water_vapor_files):
        if not i.name.startswith('.DS_Store'):
            water_vapor_data[f'wv_{i.name[-6:len(i.name)]}_data'] = read_file(i)

    #print(water_vapor_data.keys())
    #Pulling in config file for pD6 samples
    config_df = pd.read_excel(path_to_excel, sheet_name)
    #Cutting names in file_name column to match the imported files
    config_df["file_name"] = config_df["file_name"].apply(lambda file_name: f'{run_code}_{file_name[-5:len(file_name)]}') 


    #Indexing the config dataframe by file_name for simultaneous parsing with the pD6_raw_data dict below
    indexed_config_df = config_df.set_index('file_name')
    #print(indexed_pD6_config_df)

    #Initializing dict of post water vapor subtraction spectra
    cut_range_sub_wv_data = dict()


    for idx, row in indexed_config_df.iterrows():  
        print(idx)
        if idx in raw_data:
            raw_data_i = raw_data[idx]
            cut_range_sub_wv_data[f'{idx}_cut_range_sub_wv'] = cut_range_subtraction_multiple_wv(raw_data_i, water_vapor_data, row["range_start"], row["range_end"], SG_poly = 3, SG_points = 21)

    sample_raw = cut_range_sub_wv_data[sample_raw_name]
    sample_raw[0][0].plot()

    #Creating Empty Dict for second derivative of cut and subtracted data
    second_deriv_data = dict()

    #Filling it with second derivatives of all the data
    for i in cut_range_sub_wv_data:
        cut_range_sub_wv_data_i = cut_range_sub_wv_data[i]
        print(i)
        second_deriv_data[f'{i}_second_deriv'] = second_deriv(cut_range_sub_wv_data_i, show_plots=False)

    sample_second_deriv = second_deriv_data[sample_second_deriv_name]
    #print(sample_second_deriv)

    return sample_second_deriv, sample_raw
