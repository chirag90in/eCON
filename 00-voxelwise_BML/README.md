# Voxelwise BML

[00-finer_ROIs.ipynb](https://github.com/chirag90in/eCON/blob/master/00-voxelwise_BML/00-finer_ROIs.ipynb): Divides an ROI into smaller regions and outputs a new ROI mask. __righ\_insula\_10ROIs.nii.gz__ is one such output. 

[01-extract_betas.ipynb](https://github.com/chirag90in/eCON/blob/master/00-voxelwise_BML/01-extract_betas.ipynb): Extracts betas from each subject's beta map (the outputs of 3dDeconvolve), added covariate columns, and creates a dataframe for BML. Dataframe is saved a plain tab separated file in the data/ folder.

[02a-BML.r](https://github.com/chirag90in/eCON/blob/master/00-voxelwise_BML/02a-BML.r): Defines the variables of BML, runs BML, and saves the posteriors (as csv files) and the entire R workspce (file ending in .RData)

[02b-BML.stan](https://github.com/chirag90in/eCON/blob/master/00-voxelwise_BML/02b-BML.stan): Defines the BML model used in [02a-BML.r](https://github.com/chirag90in/eCON/blob/master/00-voxelwise_BML/02a-BML.r).

[03-rendering.ipynb](https://github.com/chirag90in/eCON/blob/master/00-voxelwise_BML/03-rendering.ipynb): Computes P+ from the posteriors and renders them an a brain template.
