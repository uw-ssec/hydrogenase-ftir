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
from hydrogenase_processing.second_deriv import second_deriv, first_deriv
from hydrogenase_processing.anchor_points import interact
from hydrogenase_processing.anchor_points import get_peaks, get_start_end_anchorpoints, get_all_anchor_points
from hydrogenase_processing.vaporfit import atm_subtraction
from hydrogenase_processing.vaporfit import AtmFitParams
from scipy.signal import find_peaks, peak_widths
from hydrogenase_processing.anchor_points import get_peaks, get_start_end_anchorpoints, get_all_anchor_points, baseline_spline, get_peaks_absorbance, baseline_correction#, get_peak_baseline_absorbance, plot_baseline_data
from hydrogenase_processing.peak_fit import gaussian, peak_fit, lorentzian

import ipywidgets as widgets

#addition from Eric interactive widgets
def interact(x, example_cut_sub, threshold_guess, adj_guess):
    
    #Preset threshold widget range and step by us
    threshold_widget = widgets.BoundedFloatText(
    value=threshold_guess,
    min=0,
    max=1,
    step=0.01,
    description='Threshold for peak selection(0 to 1 in 0.01 steps):',
    disabled=False
    )
    

    #Preset adj widget range and step by us
    adj_widget = widgets.BoundedFloatText(
        value=adj_guess,
        min=0,
        max=5,
        step=0.01,
        description='adj for anchor point selection(0 to 5 in 0.01 steps):',
        disabled=False
    )

    #a version of the get_all_anchor_points function that is easily inserted for the interactive widget function
    #a version of the get_peaks function that is format wise easily inserted for the interactive widget function
    def interact_with_get_peaks_and_get_all_anchor_points(threshold, adj):
        
        output = get_peaks(x, threshold) #, showplot=True)
        
        #re-extract values
        peaks_index = output[0]
        deriv_x_peak_val = output[1]
        d2ydx2_peak_val = output[2]

        wv_startIdx, wv_endIdx = get_start_end_anchorpoints(peaks_index[0], x)
        y_corr_abs = example_cut_sub[0][0].sub_spectrum
        anchor_points_raw_data = example_cut_sub[0][0].wavenb

        anchor_point_dict_output, peak_wavenumber, peak_absorbance = get_all_anchor_points(wv_startIdx, wv_endIdx, deriv_x_peak_val, anchor_points_raw_data, y_corr_abs, adj)

        #second derivative plot
        plt.subplot(2,1,1)
        plt.plot(deriv_x_peak_val, d2ydx2_peak_val, "ro",label = "peak finder peaks")
        plt.plot(x[2], x[1], label = "spline results")
        plt.title("2nd derivative plot with peak selection")
        plt.legend()

        plt.subplot(2,1,2)
        plt.plot(anchor_points_raw_data, y_corr_abs)
        plt.plot(peak_wavenumber, peak_absorbance,'ro', label='peaks')
        plt.plot(anchor_point_dict_output['wavenumber'], anchor_point_dict_output['absorbance'], 'bx', label = 'anchor_points')
        
        plt.xlabel("wavenumber")
        plt.ylabel("Absorbance")
        plt.title("")
        plt.legend()
        
        return output, anchor_point_dict_output, deriv_x_peak_val, anchor_points_raw_data, y_corr_abs
    
    #use one output because the output has to follow structure of ipywidget output and only interactive and produce non package specific objects
    interactive_results = widgets.interactive(interact_with_get_peaks_and_get_all_anchor_points, threshold = threshold_widget, adj = adj_widget)
    
    print(interactive_results)

    #break down the results
    output = interactive_results.result[0]
    anchor_point_dict_output = interactive_results.result[1]
    deriv_x_peak_val = interactive_results.result[2]
    anchor_points_raw_data = interactive_results.result[3]
    y_corr_abs = interactive_results.result[4]

    #show the output so that it's interactive
    display(interactive_results)
    
    anchor_point_dict = anchor_point_dict_output
    return anchor_point_dict, deriv_x_peak_val, anchor_points_raw_data, y_corr_abs




