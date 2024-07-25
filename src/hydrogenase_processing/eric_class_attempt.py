    #variables needed to be stored: 1. file name(as a list) 2. final peaks(dictionary? list of list?) 3. final anchor points (corresponding structure as the final peaks) 
    def interact(self, file_dir: str="../../data/opus_files/pD6", wv_dir: str="../../data/opus_files/water_vapor", threshold_init: float =0.4, adj_init: float=1.8):
        anchor_point_dict, deriv_x_peak_val, anchor_points_raw_data, y_corr_abs, file_save, threshold_save, adj_save = interact(file_dir, wv_dir, threshold_init, adj_init)
        self.saved_file = file_save
        self.anchor_point_output = anchor_point_dict
        self.thresholded_second_deriv = [deriv_x_peak_val, y_corr_abs]