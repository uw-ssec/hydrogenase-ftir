from hydrogenase_processing.cut_range import cut_range_subtraction_multiple_wv
from hydrogenase_processing.second_deriv import second_deriv
class ProSpecPy:
    def __init__(self, output_folder_path) -> None:
        self.output_folder = output_folder_path
        self.raw_data = None
        self.cut_subtracted_data = None
        self.cut_atmfitparams_obj = None
        self.cut_atmfitparameters = []
        self.batch_id = None
        self.sample_name = None
    
    def set_raw_data(self, raw_data, sample_name = None, batch_id= None):
        self.raw_data = raw_data
        self.sample_name = sample_name
        self.batch_id = batch_id

    def get_raw_data(self):
        return self.raw_data
            
    def get_subtracted_spectra(self):
        return self.cut_subtracted_data
    
    def get_atmfitparameters(self):
        return self.cut_atmfitparameters
    
    def get_atmfitparam_obj(self):
        return self.cut_atmfitparams_obj
    
    def plot_subtracted_spectra(self):
        self.get_atmfitparam_obj()[0].plot(self.sample_name, self.batch_id)
    
    def cut_range_subtract(self, raw_wv, range_start: int = 3997, range_end: int = 499, SG_poly: int = 3, SG_points: int = 21, showplot = True):
        cut_sub_data = cut_range_subtraction_multiple_wv(self.get_raw_data(), raw_wv, range_start, range_end, SG_poly, SG_points)
        # Storing the results in the object's attributes
        self.cut_atmfitparams_obj = cut_sub_data[0]
        atmfitparams_object = self.cut_atmfitparams_obj
        self.cut_atmfitparameters = atmfitparams_object[0].fit_atm_params
        self.cut_subtracted_data = cut_sub_data[1]
        if showplot:
            self.plot_subtracted_spectra()
    
    def second_derivative(self, showplots = False, save = True):
        spline_object, second_deriv_absorbance, second_deriv_wavenumber = second_deriv(self.get_atmfitparam_obj(), show_plots=showplots, sample_name=self.sample_name, batch_id=self.batch_id)
        self.second_deriv_dict = {}
        self.second_deriv_dict['UniSpline_Object'] = spline_object
        self.second_deriv_dict['absorbance'] = second_deriv_absorbance
        self.second_deriv_dict['wavenumber'] = second_deriv_wavenumber
    
    def get_second_deriv_dict(self):
        return self.second_deriv_dict
    


    

    
    

    


        
        
    