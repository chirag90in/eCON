setwd("~/Pessoa_Lab/eCON/00-voxelwise_BML")
# libraries
library(rstan)
library(parallel)

# Allocate all available cores
print(getOption("mc.cores"))
options(mc.cores = parallel::detectCores())
print(getOption("mc.cores"))

# read in data
dataset = read.csv("data/uncon_v_con_right_insula_10ROIs.txt", header = TRUE,sep = '\t')
# Convert VOX, ROI and Subj columns into factors
dataset$VOX <- factor(dataset$VOX)
dataset$ROI <- factor(dataset$ROI)
dataset$Subj <- factor(dataset$Subj)

# Print number of levels in each grouping variable
print(paste("Number of Subjects:-",nlevels(dataset$Subj)))
print(paste("Number of ROIs:-",nlevels(dataset$ROI)))
print(paste("Number of Voxels:-",nlevels(dataset$VOX)))
# Print number of cores allocated
print(paste0('Number of cores available: ', detectCores(all.tests = FALSE, logical = TRUE)))

## processing 
# Number of subjects
num_subjs <- nlevels(dataset$Subj)
# index for each subject
subj_to_idx <- 1 : num_subjs 
# Name each index with the subject's id from the dataset
names(subj_to_idx) <- levels(dataset$Subj)
idx_to_subj <- levels(dataset$Subj)
names(idx_to_subj) <- 1 : num_subjs

# Number of ROIs
num_rois <- nlevels(dataset$ROI) 
# index for each ROI
roi_to_idx <- 1 : num_rois
# Name each index with the ROI name from the dataset
names(roi_to_idx) <- levels(dataset$ROI)
idx_to_roi <- levels(dataset$ROI)
names(idx_to_roi) <- 1 : num_rois

# Number of VOX
num_vox <- nlevels(dataset$VOX)
# index for each VOX
vox_to_idx <- 1 : num_vox
# Name each index with the VOX name from the dataset
names(vox_to_idx) <- levels(dataset$VOX)
idx_to_vox <- levels(dataset$VOX)
names(idx_to_vox) <- 1 : num_vox

# format data for Stan: list if inputs that will go into the model
model_data <- list(N = nrow(dataset),
                   Y = dataset$Y,
                   bpd = dataset$BPdiff_stdz,
                   N_SUB = num_subjs,
                   N_ROI = num_rois,
                   N_VOX = num_vox,
                   sid = subj_to_idx[dataset$Subj],
                   rid = roi_to_idx[dataset$ROI],
                   vid = vox_to_idx[dataset$VOX])
# set Stan model. The model is defined in a separate stan file
model <- stan_model(file = "02b-BML.stan")
##########################################################################################
# Run the model
time0 = system.time(model_fit <- sampling(object = model, data = model_data, iter = 20000, chains = 4, cores = 4))
# Print the run time
print(time0)
# Extract posterior for the fixed (POP) and varying effects intercepts
ps_POP <- extract(model_fit,pars="POP")$POP
ps_ROI <- extract(model_fit,pars="ROI")$ROI
ps_VOX <- extract(model_fit,pars="VOX")$VOX
# Add POP posterior to VOX posteriors. 
ps_POP_VOX <- apply(ps_VOX,2,'+',ps_POP)

# Create a result directory
if (!dir.exists("03-results")){dir.create("03-results")}
# Save the POP+VOX posteriors in a csv file
ps_POP_VOX <- data.frame(ps_POP_VOX)
colnames(ps_POP_VOX) <- levels(dataset$VOX)
write.csv(ps_POP_VOX,'03-results/ps_POP_VOX.csv')

# The ROI posterior also needs to be added to POP+VOX, 
# but it is to be done more carefully. So that is done in python. 
# For now we will just save the ROI posteriors in a csv file. 
ps_ROI <- data.frame(ps_ROI)
colnames(ps_ROI) <- levels(dataset$ROI)
write.csv(ps_ROI,'03-results/ps_ROI.csv')

# Save the entire R workspace.  
save.image(file="03-results/results.RData")

