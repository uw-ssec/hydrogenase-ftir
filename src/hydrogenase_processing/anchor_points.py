
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#find peaks
from scipy.signal import find_peaks, peak_widths
from scipy.interpolate import UnivariateSpline

def get_peaks(d2ydx2_spl_upsidedown, second_deriv, threshold = 0.15, showplot = False):
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
    
    if showplot:
        plt.plot(deriv_x_peak_val, d2ydx2_peak_val, "ro",label = "peak finder peaks")
        plt.plot(second_deriv[2], second_deriv[1], label = "spline results")
        plt.legend()

    return peaks_index, deriv_x_peak_val




def get_start_end_anchorpoints(peaks_index, d2ydx2_spl_upsidedown, second_deriv):
    peak_wid = peak_widths(d2ydx2_spl_upsidedown, peaks_index, rel_height=1) 

    width_endIdx = [int(x) for x in peak_wid[2]]
    wv_endIdx =[]

    for i in width_endIdx:
        wv_endIdx.append(second_deriv[2][i])

    width_startIdx = [int(x) for x in peak_wid[3]]
    wv_startIdx = []

    for i in width_startIdx:
        wv_startIdx.append(second_deriv[2][i])
    return wv_startIdx, wv_endIdx, width_startIdx, width_endIdx


def get_all_anchor_points(wv_startIdx, wv_endIdx, deriv_x_peak_val, anchor_points_raw_data, y_corr_abs, plot_title=None, adj_factor = 1, show_plot = True):

    smaller_peak_wid = get_smaller_peak_width(deriv_x_peak_val, wv_startIdx, wv_endIdx)
    post_process_anchor_points = []
    post_process_anchor_points_abs = []

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

    if show_plot:
        #get all peak wavenumber and absorbance for plotting
        peak_wavenumber, peak_absorbance = get_peaks_absorbance(deriv_x_peak_val, anchor_points_raw_data, y_corr_abs)
        plt.plot(anchor_points_raw_data, y_corr_abs)
        plt.plot(peak_wavenumber, peak_absorbance,'ro', label='peaks')
        plt.plot(post_process_anchor_data['wavenumber'], post_process_anchor_data['absorbance'], 'bx', label = 'anchor_points')
        
        plt.title(plot_title)
        plt.xlabel("wavenumber")
        plt.ylabel("Absorbance")
        plt.legend()
        plt.plot
    return anchor_data_sorted


def get_peaks_absorbance(deriv_x_peak_val,x_wavenb, y_corr_abs):
    # Define the range
    range_width = 2
    peak_wavenumber = []
    peak_absorbance = []
    # Get the wavelengths within range of the anchor points
    wavelengths_within_range = []
    for peak_val in deriv_x_peak_val:
        indices_within_threshold = [index for index, value in enumerate(x_wavenb) if abs(value - peak_val) <= 2]
        data = pd.DataFrame({'wv': x_wavenb[indices_within_threshold], 'abs': y_corr_abs[indices_within_threshold]}) 
        peak_data = data.loc[data['abs'].idxmax()]
        peak_wavenumber.append(peak_data['wv'])
        peak_absorbance.append(peak_data['abs'])
    return peak_wavenumber, peak_absorbance


def get_smaller_peak_width(deriv_x_peak_val, wv_startIdx, wv_endIdx):
    smaller_peak_wid = []
    for i in range(len(deriv_x_peak_val)):
        left_wid = deriv_x_peak_val[i] - wv_startIdx[i]
        right_wid = wv_endIdx[i] - deriv_x_peak_val[i]
        #smaller peak width
        smaller_peak_wid.append(min(left_wid, right_wid))
    return smaller_peak_wid


def baseline_spline(anchor_points, degree=3, smooth=0):
    spline_fit = UnivariateSpline(anchor_points['wavenumber'], anchor_points['absorbance'],k = degree, s=smooth)
    x_range = np.linspace(int(min(anchor_points['wavenumber'])), int(max(anchor_points['wavenumber'])), 1000)
    baseline_fit = spline_fit(x_range)
    baseline_curve = pd.DataFrame({'wavenumber':x_range, 'absorbance': baseline_fit})
    return baseline_curve



