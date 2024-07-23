from hydrogenase_processing.cut_range import cut_range_subtraction_multiple_wv
from hydrogenase_processing.second_deriv import second_deriv
from hydrogenase_processing.anchor_points import get_peaks,get_start_end_anchorpoints,get_all_anchor_points, get_peaks_absorbance
from hydrogenase_processing.baseline import baseline_spline, baseline_correction,get_baseline_peak_index
import os
import matplotlib.pyplot as plt
import pandas as pd

class ProSpecPy:
    def __init__(self, output_folder_path= None) -> None:
        self.output_folder = output_folder_path
        if self.output_folder is not None:
            if os.path.exists(self.output_folder):
                raise FileExistsError(f"The folder '{self.output_folder}' already exists.")
            else:
                os.makedirs(self.output_folder)
        else:
            print("No output folder specified! Results from analysis will not be saved")
        self.raw_data = None
        self.cut_subtracted_data = None
        self.cut_atmfitparams_obj = None
        self.cut_atmfitparameters = []
        self.batch_id = None
        self.sample_name = None
        self.second_deriv_dict = {}
        self.second_deriv_peak_dict = {}
        self.anchor_points = None
        self.anchor_points_peak_dict = {}
        self.baseline_curve = None
    
    def set_raw_data(self, raw_data, sample_name = None, batch_id= None):
        self.raw_data = raw_data
        self.sample_name = sample_name
        self.batch_id = batch_id

    def get_raw_data(self):
        return self.raw_data
    
    def save_plot(self, fig, filename, verbose=True):
        """
        Save the given plot to the specified output folder.

        Args:
            fig (matplotlib.figure.Figure): The figure object to save.
            filename (str): The name of the file to save the figure as.
        """
        if not os.path.exists(self.output_folder):
            print("No such folder exists! please create a valid output folder")
        elif self.output_folder == None:
            print("No output folder given")
        else:
            full_path = os.path.join(self.output_folder, filename)
            fig.savefig(full_path)
            if verbose:
                print(f"Plot saved to {full_path}")
            
    def get_subtracted_spectra(self):
        return self.cut_subtracted_data
    
    def get_atmfitparameters(self):
        return self.cut_atmfitparameters
    
    def get_atmfitparam_obj(self):
        return self.cut_atmfitparams_obj
    
    def get_subtracted_spectra_absorbance(self):
        return self.get_atmfitparam_obj()[0].sub_spectrum
    
    def get_subtracted_spectra_wavenumber(self):
        return self.get_atmfitparam_obj()[0].wavenb
    
    def plot_subtracted_spectra(self, save = True, showplots = False):
        subtracted_fig = self.get_atmfitparam_obj()[0].plot(self.sample_name, self.batch_id, showplots =showplots)
        return subtracted_fig
    
    def cut_range_subtract(self, raw_wv, range_start: int = 3997, range_end: int = 499, SG_poly: int = 3, SG_points: int = 21, showplot = False, save = True, verbose=True):
        cut_sub_data = cut_range_subtraction_multiple_wv(self.get_raw_data(), raw_wv, range_start, range_end, SG_poly, SG_points)
        # Storing the results in the object's attributes
        self.cut_atmfitparams_obj = cut_sub_data[0]
        atmfitparams_object = self.cut_atmfitparams_obj
        self.cut_atmfitparameters = atmfitparams_object[0].fit_atm_params
        self.cut_subtracted_data = cut_sub_data[1]
        cut_subtract_fig = self.plot_subtracted_spectra(showplots = showplot)
        if save:
            filename = 'cut_range_subtracted_spectra'
            self.save_plot(cut_subtract_fig,filename, verbose=verbose)
            # Save cut_atmfitparameters as CSV
            atmfitparams_df = pd.DataFrame([self.cut_atmfitparameters])
            atmfitparams_filename = 'cut_atmfitparameters.csv'
            atmfitparams_df.to_csv(os.path.join(self.output_folder, atmfitparams_filename), index=False)
            if verbose:
                print(f"Atmospheric fit parameters saved to {os.path.join(self.output_folder, atmfitparams_filename)}")
            
            # Save cut_subtracted_data as CSV
            cut_subtracted_data_df = pd.DataFrame(self.cut_subtracted_data)
            cut_subtracted_data_filename = 'cut_subtracted_data.csv'
            cut_subtracted_data_df.to_csv(os.path.join(self.output_folder, cut_subtracted_data_filename), index=False)
            if verbose:
                print(f"Cut subtracted data saved to {os.path.join(self.output_folder, cut_subtracted_data_filename)}")
    
    def second_derivative(self, showplots = False, save = True, verbose = True):
        self.second_deriv_tuple, cut_subtracted_data_fig, second_derivative_fig = second_deriv(self.get_atmfitparam_obj(), show_plots=showplots, sample_name=self.sample_name, batch_id=self.batch_id)
        self.second_deriv_dict['UniSpline_Object'] = self.second_deriv_tuple[0]
        self.second_deriv_dict['absorbance'] = self.second_deriv_tuple[1]
        self.second_deriv_dict['wavenumber'] = self.second_deriv_tuple[2]
        if save:
            filename = 'subtracted_spectra'
            self.save_plot(cut_subtracted_data_fig,filename, verbose=verbose)
            if verbose:
                print(f"subtracted_spectra plot saved to {os.path.join(self.output_folder, filename)}")

            second_deriv_filename = 'second_derivative_fig'
            self.save_plot(second_derivative_fig,second_deriv_filename, verbose=verbose)
            if verbose:
                print(f"second derivative plot saved to {os.path.join(self.output_folder, second_deriv_filename)}")

            data_df = pd.DataFrame({
                'wavenumber': self.second_deriv_dict['wavenumber'],
                'absorbance': self.second_deriv_dict['absorbance']
            })
            csv_filename = 'second_derivative_data.csv'
            data_df.to_csv(os.path.join(self.output_folder, csv_filename), index=False)
            if verbose:
                print(f"Second derivative csv data saved to {os.path.join(self.output_folder, csv_filename)}")

    def get_second_deriv_dict(self):
        return self.second_deriv_dict
    
    def get_second_deriv_tuple(self):
        return self.second_deriv_tuple

    def peak_finder(self, threshold):
        self.threshold = threshold
        self.peak_information, peak_wavenumber, peak_seconderiv_absorbance = get_peaks(self.get_second_deriv_tuple(), self.threshold)
        peak_index = self.peak_information[0]
        self.second_deriv_peak_dict['peak_index'] = peak_index
        self.second_deriv_peak_dict['peak_wavenumber'] = peak_wavenumber
        self.second_deriv_peak_dict['peak_second_deriv_absorbance'] = peak_seconderiv_absorbance

    def get_second_deriv_peak_dict(self):
        return self.second_deriv_peak_dict
    
    def get_peak_index(self):
        return self.second_deriv_peak_dict['peak_index']
    
    def save_second_deriv_peak_plot(self):
        #TO DO
        pass
        
    def anchor_point_fit(self, adjustment_factor):
        wv_startIdx, wv_endIdx = get_start_end_anchorpoints(self.get_peak_index(), self.get_second_deriv_tuple())
        anchor_points, peak_wavenumber,peak_absorbance = get_all_anchor_points(wv_startIdx, wv_endIdx, self.second_deriv_peak_dict['peak_wavenumber'], self.get_subtracted_spectra_wavenumber(), self.get_subtracted_spectra_absorbance(), adjustment_factor)
        self.anchor_points = anchor_points
        self.anchor_points_peak_dict['peak_wavenumber'] = peak_wavenumber
        self.anchor_points_peak_dict['peak_absorbance'] = peak_absorbance
    
    def get_anchor_points_peak_dict(self):
        return self.anchor_points_peak_dict
    
    def get_anchor_points(self):
        return self.anchor_points
    
    def baseline_fit(self):
        self.baseline_curve = baseline_spline(self.get_anchor_points())
    
    def get_baseline_curve(self):
        return self.baseline_curve
    
    def subtract_baseline(self):
        self.baseline_corrected_abs = baseline_correction(self.get_baseline_curve(), self.get_subtracted_spectra_wavenumber(),self.get_subtracted_spectra_absorbance())
        peak_wv, peak_abs = get_peaks_absorbance(self.second_deriv_peak_dict['peak_wavenumber'], self.get_subtracted_spectra_wavenumber(),self.get_subtracted_spectra_absorbance())
        peak_wv_index, peak_wv_baseline, peak_baseline_abs = get_baseline_peak_index(self.baseline_corrected_abs, self.get_subtracted_spectra_wavenumber(), peak_wv) 
        self.baseline_corrected_peak_dict['peak_index'] = peak_wv_index
        self.baseline_corrected_peak_dict['wavenumber'] = peak_wv_baseline
        self.baseline_corrected_peak_dict['absorbance'] = peak_baseline_abs

        

    
    

    


        
        
    