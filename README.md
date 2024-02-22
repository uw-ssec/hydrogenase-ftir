# Hydrogenase FTIR Data Processing #



Pure Python code attempting to replicate and package the following data processing steps:


The manual process in OPUS:  
- Subtract water vapour (subtract a “reference” water vapour spectra) 			 
- Cut data range to area of interest 										 

These files are then imported into Origin as .dpt files to: 						 
- Plot data 	 
- Take second derivative (of Y-variable, absorbance) and plot  
- Add individual points to define baseline – an iterative and subjective process relying on points identified from the second derivative plot 					 
- Fit peaks from your baseline-subtracted plot

## Current Status: 12/12/23 ##
First Step Replicated in [`opus_import.ipynb`](src/opus_import.ipynb), using subtraction algorithm found in [`vaporfit.ipynb`](src/vaporfit.ipynb)

Data Files can be found in [`opus_files`](./opus_files) folder, with plots in the [`output_plots`](./output_plots) folder

