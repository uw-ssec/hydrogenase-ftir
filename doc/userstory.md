## User story ##
1. Elizebeth is a research scientist. Elizebeth wants to […]. She uses the software Origin and OPUS often. Elizebeth wants a quick and easy software for peak analysis. Elizebeth’s job involves analyzing the spectra data and making sure they are scientifically sound, she values the accuracy and speed of the spectra analysis.

Note: […] are any use cases written below


2. Reachel is a student in Elizebeth’s lab. Reachel wants to […]. She had used the Origin and OPUS software a few times. Elizebeth wants a quick and easy software for peak analysis. Elizebeth’s job involves analyzing the spectra data and making sure they are scientifically sound, she values the accuracy, speed of the spectra analysis, and how easy the software is to be picked up by someone new.

Note: […] are any use cases written below

3. Joe is a journal reviewer. Joe wants to understand what has been done in the paper. He had never used Origin or OPUS. Joe wants easy to understand model for spectra analysis. Joe’s jobs involves understanding software and their scientific soundness, he values accuracy, simplicity, and efficiency.






## Use Cases (according to the video description of our researcher collaborator) ##

Subtract water vapour (subtract a “reference” water vapour spectra) – 
- 3:21-4:50
Steps: 
1.	Open water vapor file form folder
2.	Import water vapor file into OPUS
3.	 Select AB spectrum option under the file, click on a display option
4.	Open data file form folder
5.	Import data file into OPUS individually to avoid random ordering of the files such as 0,1,2,3,4 gets arranged as 4,2,3,1,0 after imported.
- 7:20-9:25
1. Open the last spectrum of a collection of volage potential (naming convention: number+a/b/c(different potentials)+ as iso Hyd2 dark titration+ mV(potential)+.number of run) (spectrum of same potentials are used to check whether the experiment had reached equilibrium)
2. Zoom into 1800-2200, interested species in 2100-2050, 1980-1880
3. Difficult to tell the spectra apart when they are yellow or green lines
4. Shift curve-> whole curve to put the data curves in order of the potential, where most positive potential is at the top 


- 9:25-11:34
Actual Subtraction Step 
Steps:
1.	Zoom in to where water vapor is most seen within the 1800-2200 range
2.	Water vapor times a 1-100 value and choose which digit to subtract eg 0.1 or 0.001
3.	Store after removing water vapor to satisfaction. Appropriate peak seems to involve not introducing positive peaks in the water vapor 10:45 (final results can be changed later on and goal seems to remove secondary noise peaks )
4.	Save file as original name + _subWV.dpt (save file No.1)



Cut data range to area of interest 11:34-13:00	
1.	Cut to 1800 to 2150 (desert the rest of the data)
2.	Create a new save file named original name+_subWV cut.dpt (save file No.2)									 
These files are then imported into Origin as .dpt files: 13:00-15:42
1.	Drag the cut and uncut data into origin. A promote show up No match filters found or the D&D feature of the related filter is turned off. Use import wizard to import the 13:17 		
2.	Worksheet with file name and workbook with file name, so you don’t have to go back and rename files
3.	Put cut data and uncut data into separate folders
4.	Might ant look different frequency for D2O(15:34)
5.	Make individual folder for each potential
			  
Plot data, 17:30-23:15 
Desired plot format and looks as at 22:29


Take second derivative (of Y-variable, absorbance) and plot  23:15-25:45 
1.	Recalculate as manual To avoid origin doing long time auto update for a change in the data
2.	Take second derivate and output as a new column
3.	Plot it for different data, and name them 
4.	Display peak wave number values
5.	Analyse any peak that is higher than noise magnitude, any peak in the right speice range of the same noise magnitude also would analyse, but might want to take a grain of salt, need a label (2 interest region CO and CN)
6.	write down where all the peaks that will be analysed with their corresponding wavenumber
7.	wants optional display of wavenumber on the final plot in 6 (28:06)

Add individual points to define baseline – an iterative and subjective process relying on points identified from the second derivative plot,25:45-28:07,  			

(28:10) Define baseline:
1.	select absorbance, manual (no auto bec origin tries to do too much and leads to crashing)
2.	user defined baseline, use 2nd derivative peaks (snap to spectrum 29:29)
3.	origin glitch where smoothing on the curve leave 8 2nd derivative zeros in random spots, corrections will be needed and will need more than 8 points (always tell program to find 500 points, so that you can always add more points onto it)
4.	Redefine where the points are:
-	Connect the points by interpolation (interpolation method: spline)
-	Add points at the very end of the spectrum on both ends
-	proceed 31:43 to put 17 data points onto the curve in a roughly on the curve manner, make sure these are no where near the peaks noted down, or you will drop the peak out, marked area 1800-1840 And 2120-2160, no anker points on the peak around 50 wavenumber is okay tho for example the 2 points shown here should be removed
 
-	Look back at the 2nd derivative (used to tell where peaks are) plot, use the peak wavenumber as location to, anker points on both sides of that peak, try not to put points betw peaks bec you can manually shorten the peak tail
-	Zoom in to remove negative peak so no points below baseline
-	Then select auto subtract and auto rescale
-	Check results to match peaks of subtracted and 2nd derivative peaks, add more anker points to show peaks that were supposed to show
-	Use 4-5% threshold(x% or more of peak height gets identified as peak) to calculate threshold peak hight, do not want to find peaks that are shoulders
-	Probably wants a progress bar, they seem to be scared of a spinning loading icon lol
-	41:05 rough anker point plot
 
-	Fit the data with Gaussian and Gaussian LorenCross and lorencian (later 2 used for peaks with really long tails)
- 44:30 don’t want very narrow sharp or really wide peaks. FWHM 4-5 is about alright
-	Modification required when
 
Black line- actual data
Blue line- individual peak fit
Red line- overall cumulative peak fit
-	Correct the red line for where you think the maximum of the peak is at
-	Add a peak where missing, however if adding lowers other peaks:
o	Then add in anker point between peaks
-	If peak is really wide
o	Go back and readjust the anker point
-	Generate a report
-	Gaussian LorenCross need peak center position
Fit peaks from your baseline-subtracted data (informed by second derivative plot)28:07-54:26 

Plot subtracted baseline data  1:02:00-1:05:00 
-	Plot subtracted baseline data
-	Do a second derivative on the subtracted baseline data compare to 2nd derivative of raw data
-	Match the 2 second derivative plots, noise does not matter just main peaks need to match

## Components ##
[only some of the use cases are solved with components at the current state, others are either too complex to complete within the time span]
1.	Input data storage folder or dataframe
2.	A sorting function such that the input data are listed in order
3.	A subtraction percentage function that times the water vapor by an adjustment factor and subtract that value from the data
4.	A cutting function that cuts the subtracted data to wavenumber 1800-2150, and output a save file with original name_subWV cut.dpt as name
5.	A sorting function that put all cut and uncut savefile into separate folders (in future developement)
6.	A sorting function that put different potential data into separate folders (in future development)
7.	A function that takes second derivative of the data
8.	A function that allows users to choose what is noise and irrelevant data and what is included as signal peaks
9.	A data object that records where each peak is at
10.	An interactive function to assign anchor points on the data
11.	A function subtract the baseline established by the anchor point
12. A function that use Guassian and Lorenzian fitting for the peaks left after the baseline subtraction and present result visually

