Light Beads Microscopy with MAxiMuM Module: Data processing pipeline and work flow.

The following softwares are required to be installed and added to the MATLAB directory prior to running this pipeline:
- CaImAn (MATLAB)
- NoRMCorrE (MATLAB)
- Scanimage (MATLAB)

Work flow:

After saving data in raw .tif files, run preProcessMAxiMuM to reconstruct each frame for each plane from consitituent ROIs and save the data set as a series of separate files with corresponding motion-corrected videos of each plan.

Next, run planarSegmentation to segment each plane and find neurons. 

Run calculateOffset to refine estimates of the lateral shifts between each plane.

Finally, run comparePlanes to collate the results from each plane into a single file containing all the traces and coordinates of the extracted neurons. 

