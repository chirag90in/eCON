# Voxelwise BML

[00-finer_ROIs.ipynb](https://github.com/chirag90in/eCON/blob/master/00-voxelwise_BML/00-finer_ROIs.ipynb): Divides an ROI into smaller regions and outputs a new ROI mask. __righ\_insula\_10ROIs.nii.gz__ is one such output. 

[01-extract_betas.ipynb](https://github.com/chirag90in/eCON/blob/master/00-voxelwise_BML/01-extract_betas.ipynb): Extracts betas from each subject's beta map (the outputs of 3dDeconvolve), added covariate columns, and creates a dataframe for BML. Dataframe is saved a plain tab separated file in the data/ folder.
