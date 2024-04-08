
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#find peaks
from scipy.signal import find_peaks, peak_widths


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
    #sec_dev_endIdx = []
    for i in width_endIdx:
        wv_endIdx.append(second_deriv[2][i])
        #sec_dev_endIdx.append(x[1][i])


    width_startIdx = [int(x) for x in peak_wid[3]]
    wv_startIdx = []
    #sec_dev_startIdx = []
    for i in width_startIdx:
        wv_startIdx.append(second_deriv[2][i])
        #sec_dev_startIdx.append(x[1][i])

    '''anchor_points_from_second_deriv = []
    for wavenb in wv_startIdx:
        anchor_points_from_second_deriv.append(wavenb)
    for wavenb in wv_endIdx:
        anchor_points_from_second_deriv.append(wavenb)

    anchor_points_from_second_deriv = set(anchor_points_from_second_deriv)'''

    return wv_startIdx, wv_endIdx, width_startIdx, width_endIdx


def get_all_anchor_points(wv_startIdx, wv_endIdx, deriv_x_peak_val, anchor_points_raw_data, y_corr_abs, plot_title=None, adj_factor = 1):
    #for loop to compute the left hand side and right hand side peak width for each peak, then put them into a list
    peak_wid_2_sides = []
    for i in range(len(deriv_x_peak_val)):
        left_wid = deriv_x_peak_val[i] - wv_startIdx[i]
        right_wid = wv_endIdx[i] - deriv_x_peak_val[i]
        peak_wid_2_sides.append([left_wid, right_wid])

    #print(peak_wid_2_sides)

    #for loops to decide which anchor point should be included becuase it is far enought away from the peak
    post_process_anchor_points = []
    post_process_anchor_points_abs = []

    peak_correspondence = []
    anchorpt_peak_correspondence = []

    # Peaking sorting: subtract each anchor point wi the peak and compare which difference is the smallest among the peaks and this one anchor point to figure out which peak this anchor point is close to
    for index in range(len(anchor_points_raw_data)):
        #print(anchor_points_raw_data[index])
        peak_to_anchor = []
        # for every peak
        for i in range(len(deriv_x_peak_val)):
            #print('peak',i, 'is', deriv_x_peak_val[i])
            #peak sorting
            pa_distance = abs(deriv_x_peak_val[i] - anchor_points_raw_data[index])
            peak_to_anchor.append(pa_distance)
        index_of_min = peak_to_anchor.index(min(peak_to_anchor))
        #print(anchor_points_raw_data[index], 'is closer to peak', index_of_min)
        peak_correspondence.append(index_of_min)
    anchorpt_peak_correspondence = [anchor_points_raw_data, peak_correspondence]

            
    #for loop: 1st: determine is the anchor point on the left or right of the corresponding peak, 2nd compare the difference between the anchor points and the peak against the peak width stored in peak_wid_2_sides to determine if they are far enough away from the peak         
    for index in range(len(anchor_points_raw_data)):
            corr_peak_index = anchorpt_peak_correspondence[1]
            if anchor_points_raw_data[index] < deriv_x_peak_val[corr_peak_index[index]]:
                #print('on the left of', deriv_x_peak_val[anchorpt_peak_correspondence[1][index]])
                if deriv_x_peak_val[corr_peak_index[index]] - anchor_points_raw_data[index] > peak_wid_2_sides[corr_peak_index[index]][0]*adj_factor:
                    #print('include')
                    post_process_anchor_points.append(anchor_points_raw_data[index])
                    post_process_anchor_points_abs.append(y_corr_abs[index])
            else:
                #print('on the right of', deriv_x_peak_val[i])
                if anchor_points_raw_data[index] - deriv_x_peak_val[corr_peak_index[index]] > peak_wid_2_sides[corr_peak_index[index]][1]*adj_factor:
                    #print('include')
                    post_process_anchor_points.append(anchor_points_raw_data[index])
                    post_process_anchor_points_abs.append(y_corr_abs[index])


    #post processesing to avoid repeating values and make sure the wavenumbers are in the same acending or decending order
    post_process_anchor_data = pd.DataFrame({'wavenumber': post_process_anchor_points, 'absorbance': post_process_anchor_points_abs})
    post_process_anchor_data = post_process_anchor_data.drop_duplicates()
    anchor_data_sorted = post_process_anchor_data.sort_values(by='wavenumber').reset_index()
    #print(anchor_data_sorted)

    #get all peak wavenumber and absorbance for plotting
    peak_wavenumber, peak_absorbance = get_peaks_absorbance(deriv_x_peak_val, anchor_points_raw_data, y_corr_abs)
    plt.plot(anchor_points_raw_data, y_corr_abs)
    plt.plot(peak_wavenumber, peak_absorbance,'ro')
    plt.plot(post_process_anchor_data['wavenumber'], post_process_anchor_data['absorbance'], 'bx')
    
    plt.title(plot_title)
    plt.xlabel("wavenumber")
    plt.ylabel("Absorbance")
    plt.plot


def get_peaks_absorbance(deriv_x_peak_val,x_wavenb, y_corr_abs):
    # Define the range
    range_width = 2
    peak_wavenumber = []
    peak_absorbance = []
    # Get the wavelengths within range of the anchor points
    wavelengths_within_range = []
    for peak_val in deriv_x_peak_val:
        indices_within_threshold = [index for index, value in enumerate(x_wavenb) if abs(value - peak_val) <= 2]
        for idx in indices_within_threshold:
            #peak_wv - one_end = peak_wid
            #if x_wavenb[idx] <= peak_wid[0] + peak_wv:
            peak_wavenumber.append(x_wavenb[idx])
            peak_absorbance.append(y_corr_abs[idx])
    return peak_wavenumber, peak_absorbance