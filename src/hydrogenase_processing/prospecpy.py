from hydrogenase_processing.cut_range import cut_range_subtraction_multiple_wv
import numpy as np
class ProSpecPy:
    def __init__(self, output_folder_path) -> None:
        self.output_folder = output_folder_path
        self.raw_data = None
        self.cut_subtracted_data = None
        self.cut_atmfitparams_obj = None
    
    def set_raw_data(self, raw_data):
        self.raw_data = raw_data

    def get_raw_data(self):
        return self.raw_data
    
    def cut_range_subtract(self, raw_wv, range_start: int = 3997, range_end: int = 499, SG_poly: int = 3, SG_points: int = 21):
        """
        Cuts the specified range from raw spectra data and performs atmospheric subtraction using multiple water vapor spectra.

        Args:
            raw_wv (dict): Dict of all water vapor spectra, with appropriately named keys.
            range_start (int, optional): Start wavenumber for the desired range. Defaults to 3997.
            range_end (int, optional): End wavenumber for the desired range. Defaults to 499.
            SG_poly (int, optional): Polynomial order for Savitzky-Golay smoothing. Defaults to 3.
            SG_points (int, optional): Number of points for Savitzky-Golay smoothing. Defaults to 21.

        Returns:
            Tuple[List['AtmFitParams'], np.ndarray]: Atmospheric subtraction results.
                - List of atmospheric fitting parameters (e.g., coefficients for water vapor correction).
                - NumPy array representing the results of atmospheric subtraction.
        """
        cut_sub_data = cut_range_subtraction_multiple_wv(self.get_raw_data(), raw_wv, range_start, range_end, SG_poly, SG_points)
        
        # Storing the results in the object's attributes
        self.cut_atmfitparams_obj = cut_sub_data[0]
        self.cut_subtracted_data = cut_sub_data[1]
        
    def get_subtracted_spectra(self):
        return self.cut_subtracted_data
    
    def get_atmfitparams_object(self):
        return self.cut_atmfitparams_obj

    


        
        
    