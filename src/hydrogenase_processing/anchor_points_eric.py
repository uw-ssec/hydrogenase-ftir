
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#find peaks
from scipy.signal import find_peaks, peak_widths
from scipy.interpolate import UnivariateSpline
#widget
import ipywidgets as widgets


def get_peaks(second_deriv, threshold = 0.15): #, showplot = False):
    """
    Function to detect peaks in the second derivative of a spline function.

    Parameters:
    - second_deriv: tuple
        Tuple containing the results of the spline function, including x and y values.
    - threshold: float, optional (default=0.15)
        Threshold value used to determine the prominence of peaks.
    - showplot: bool, optional (default=False)
        Flag to indicate whether to plot the detected peaks.

    Returns:
    - peaks_index: array-like
        Indices of the detected peaks.
    - deriv_x_peak_val: array-like
        x-coordinate values(wavenumbers) of the detected peaks.
    """
    d2ydx2_spl_upsidedown = second_deriv[1] * -1
    relative_height = threshold * max(d2ydx2_spl_upsidedown)
    peaks_index = find_peaks(d2ydx2_spl_upsidedown, prominence=relative_height)

    #use for loops to extract the coordinates of the peaks so we can plot them on the plot above
    d2ydx2_peak_val = []
    deriv_x_peak_val = []

    for i in peaks_index[0]:
        d2ydx2_peak = second_deriv[1][i]
        deriv_x_peak = second_deriv[2][i]

        d2ydx2_peak_val.append(d2ydx2_peak)
        deriv_x_peak_val.append(deriv_x_peak)

# inactivated for the plot to be in desired layout in the interact 
#    if showplot:    
#        plt.plot(deriv_x_peak_val, d2ydx2_peak_val, "ro",label = "peak finder peaks")
#        plt.plot(second_deriv[2], second_deriv[1], label = "spline results")
#        plt.legend()

#added d2ydx2_peak_val for ease of plotting data navigation
    return peaks_index, deriv_x_peak_val, d2ydx2_peak_val




def get_start_end_anchorpoints(peaks_index, second_deriv):

    """
    Function to determine the start and end anchor points of peaks.

    Parameters:
    - peaks_index: array-like
        Indices of the detected peaks.
    - second_deriv: tuple
        Second derivative of the spline function

    Returns:
    - wv_startIdx: list
        x-coordinate values(wavenumber) of the start anchor points of peaks.
    - wv_endIdx: list
        x-coordinate values(wavenumber) of the end anchor points of peaks.
    - width_startIdx: list
        Indices of the start anchor points of peaks.
    - width_endIdx: list
        Indices of the end anchor points of peaks.
    """
    d2ydx2_spl_upsidedown = second_deriv[1] * -1
    peak_wid = peak_widths(d2ydx2_spl_upsidedown, peaks_index, rel_height=1) 

    width_endIdx = [int(x) for x in peak_wid[2]]
    wv_endIdx =[]

    for i in width_endIdx:
        wv_endIdx.append(second_deriv[2][i])

    width_startIdx = [int(x) for x in peak_wid[3]]
    wv_startIdx = []

    for i in width_startIdx:
        wv_startIdx.append(second_deriv[2][i])
    return wv_startIdx, wv_endIdx

#plot title deleted
def get_all_anchor_points(wv_startIdx, wv_endIdx, deriv_x_peak_val, anchor_points_raw_data, y_corr_abs, adj_factor=1): #, show_plot = True, plot_title=None,):
    """
    Function to filter and post-process anchor points based on peak characteristics.

    Parameters:
    - wv_startIdx: list
        x-coordinate values (wavenumber) of the start anchor points of peaks.
    - wv_endIdx: list
        x-coordinate values (wavenumber) of the end anchor points of peaks.
    - deriv_x_peak_val: array-like
        x-coordinate values(wavenumber) of the detected peaks.
    - anchor_points_raw_data: array-like
        Raw spectra data.
    - y_corr_abs: array-like
        Corresponding absorbance values.
    - plot_title: str, optional (default=None)
        Title for the plot.
    - adj_factor: float, optional (default=1)
        Adjustment factor for filtering anchor points.
    - show_plot: bool, optional (default=True)
        Flag to indicate whether to display the plot.

    Returns:
    - anchor_data_sorted: DataFrame
        DataFrame containing sorted anchor points data.
    """
    #get the smaller width for each peak using get_smaller_peak_width()
    smaller_peak_wid = get_smaller_peak_width(deriv_x_peak_val, wv_startIdx, wv_endIdx)
    post_process_anchor_points = []
    post_process_anchor_points_abs = []
    """
    Algorithm:
        For each point in the raw data:
            1. get the index of the peak that is closest to this point using the criteria minimum(abs(peak_wavenumber - raw_point))
            2. if abs(closest peak - raw data point) > smaller peak width * adjustment factor:
                        then append the raw data point as anchor point and its absorbance 
    """

    for index in range(len(anchor_points_raw_data)):
        dist_peak_to_anchor = abs(deriv_x_peak_val-anchor_points_raw_data[index])
        closest_peak_idx = np.argmin(dist_peak_to_anchor)

        if abs(deriv_x_peak_val[closest_peak_idx] - anchor_points_raw_data[index]) > smaller_peak_wid[closest_peak_idx]*adj_factor:
                post_process_anchor_points.append(anchor_points_raw_data[index])
                post_process_anchor_points_abs.append(y_corr_abs[index])

    #post processesing to avoid repeating values and make sure the wavenumbers are in the same acending or decending order
    post_process_anchor_data = pd.DataFrame({'wavenumber': post_process_anchor_points, 'absorbance': post_process_anchor_points_abs})
    post_process_anchor_data = post_process_anchor_data.drop_duplicates()
    anchor_data_sorted = post_process_anchor_data.sort_values(by='wavenumber').reset_index()

#    if show_plot:
    #get all peak wavenumber and absorbance for plotting
    peak_wavenumber, peak_absorbance = get_peaks_absorbance(deriv_x_peak_val, anchor_points_raw_data, y_corr_abs)
#        plt.plot(anchor_points_raw_data, y_corr_abs)
#        plt.plot(peak_wavenumber, peak_absorbance,'ro', label='peaks')
#        plt.plot(post_process_anchor_data['wavenumber'], post_process_anchor_data['absorbance'], 'bx', label = 'anchor_points')
        
#        plt.title(plot_title)
#        plt.xlabel("wavenumber")
#        plt.ylabel("Absorbance")
#        plt.legend()
#        plt.plot
    return anchor_data_sorted, peak_wavenumber, peak_absorbance


def get_peaks_absorbance(deriv_x_peak_val,x_wavenb, y_corr_abs):

    """
    Function to retrieve peak wavenumbers and corresponding absorbance values.

    Parameters:
    - deriv_x_peak_val: array-like
        x-coordinate values (wavenumber) of the detected peaks.
    - x_wavenb: array-like
        Wavenumber values of raw spectra.
    - y_corr_abs: array-like
        Corresponding absorbance values of raw spectra.

    Returns:
    - peak_wavenumber: list
        Wavenumber values corresponding to the peaks.
    - peak_absorbance: list
        Absorbance values corresponding to the peaks.
    """
    # Define the range
    range_width = 2
    peak_wavenumber = []
    peak_absorbance = []

    for peak_val in deriv_x_peak_val:
        indices_within_threshold = [index for index, value in enumerate(x_wavenb) if abs(value - peak_val) <= range_width]
        data = pd.DataFrame({'wv': x_wavenb[indices_within_threshold], 'abs': y_corr_abs[indices_within_threshold]}) 
        peak_data = data.loc[data['abs'].idxmax()]
        peak_wavenumber.append(peak_data['wv'])
        peak_absorbance.append(peak_data['abs'])
    return peak_wavenumber, peak_absorbance


def get_smaller_peak_width(deriv_x_peak_val, wv_startIdx, wv_endIdx):
    """
    Function to calculate the smaller peak width for each peak.

    Parameters:
    - deriv_x_peak_val: array-like
        x-coordinate values (wavenumber) of the detected peaks.
    - wv_startIdx: list
        x-coordinate values (wavenumber) of the start anchor points of peaks.
    - wv_endIdx: list
        x-coordinate values (wavenumber) of the end anchor points of peaks.

    Returns:
    - smaller_peak_wid: list
        Smaller peak width for each peak.
    """
    smaller_peak_wid = []
    for i in range(len(deriv_x_peak_val)):
        left_wid = deriv_x_peak_val[i] - wv_startIdx[i]
        right_wid = wv_endIdx[i] - deriv_x_peak_val[i]
        #smaller peak width
        smaller_peak_wid.append(min(left_wid, right_wid))
    return smaller_peak_wid


def baseline_spline(anchor_points, degree=3, smooth=0):
    """
    Function to fit a spline curve to anchor points data to estimate baseline.

    Parameters:
    - anchor_points: DataFrame
        DataFrame containing anchor points with 'wavenumber' and 'absorbance' columns.
    - degree: int, optional (default=3)
        Degree of the spline interpolation.
    - smooth: float, optional (default=0)
        Smoothing parameter for spline fitting.

    Returns:
    - baseline_curve: DataFrame
        DataFrame containing the fitted baseline curve with 'wavenumber' and 'absorbance' columns.
    """
    spline_fit = UnivariateSpline(anchor_points['wavenumber'], anchor_points['absorbance'],k = degree, s=smooth)
    x_range = np.linspace(int(min(anchor_points['wavenumber'])), int(max(anchor_points['wavenumber'])), 1000)
    #x_range = anchor_points['wavenumber']
    baseline_fit = spline_fit(x_range)
    baseline_curve = pd.DataFrame({'wavenumber':x_range, 'absorbance': baseline_fit})
    return baseline_curve

def baseline_correction(baseline_points, raw_wavenumber, raw_absorbance):
    """
    Perform baseline correction on raw absorbance data using baseline points.

    Args:
    - baseline_points (DataFrame): DataFrame containing baseline wavenumber and absorbance values.
    - raw_wavenumber (array): List of wavenumber values from raw data.
    - raw_absorbance (array): List of absorbance values from raw data.

    Returns:
    - baseline_corrected_abs (list): List of baseline-corrected absorbance values.
    """
    baseline_corrected_abs = []
    for idx, wv_num in enumerate(raw_wavenumber): #iterate through each data point in raw data
        #subtract the wv_num from the baseline[wavenumber] and take absolute of the diff and sort it and get the
        #idx of the point from the baseline that is closest to the raw datapoint
        diff_array = abs(baseline_points['wavenumber'] - wv_num)
        #this is the index of the closest wavenumber in the baseline curve to that of the datapoint
        closest_wv_num = diff_array.idxmin()
        #Now subtract the baseline absornace from the raw data absorbance
        raw_minus_baseline = raw_absorbance[idx] - baseline_points.loc[closest_wv_num, 'absorbance']

        if raw_minus_baseline <0:
            #if the difference is negative, then baseline point is higher than raw absorbance which is not possible. Hence appending 0 at those points
            baseline_corrected_abs.append(0)
        else:
            baseline_corrected_abs.append(raw_minus_baseline)

            
    return baseline_corrected_abs

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


