This folder holds a functioning script for particle analysis, example images, a custom ROI for image analysis, and an example of data after workup in a csv.

The script: ratio_image_analysis_Boyt_Greytak.m is divided into sections which should be ran in sequence. Simple modifications must be made to ensure the correct data is analyzed, these are as follows:

Line 43: path to input file (folder of images)

Line 46: path to output file (location where output text file will be created)

Line 139: path to region of interest (ROI). The ROI used in our manuscripts, myroi.mat, is provided for use. This ROI is 1920x2560. Customize this size to fit your preferred ROI, or forgo an ROI entirely by skipping the "Load an ROI" section (lines 136-148)

Line 154: adjust for the amount of input files in folder, or which files are desired for analysis. With an input folder of 20 jpgs, "1:20" will analyze all files, while "1:10" will only analyze the first 10.






After the "compute statistics for each slice" section is ran, files will be analyzed individually, and a text file "particles.txt" will be generated in the output folder pathed to on line 46. This file is best opened using Excel, delimited, separated by spaces. 

All data can be sorted by "Valid", where a value of 1 denotes a particle of interest, and a value of 0 denotes a "particle" which is likely outside of the ROI, irregular in shape, irregular in luminosity, overlapping with other particles, a physical contaminant such as dust, or a spectral contaminant in the imaging device. 

This sheet provides statistics on brightness in RBG and overall brightness, standard deviation of this brightness across the individual particle, the brightest value on the particle, and the area, perimeter, and circumference of the particle in relative units. An example file was created for the sample data provided.