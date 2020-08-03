
library(brms)
library(dplyr)
library(parallel)

dataset <- read.table('data/uncon_v_con_right_insula_10ROIs.txt',header = TRUE,sep = "\t")
head(dataset)

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

# Set nuber of cores to use. 
# To run each chain on a single core, set number of core to 4
print(getOption("mc.cores"))
options(mc.cores = parallel::detectCores())
print(paste0('Number of cores allocated',getOption("mc.cores")))

# Number of iterations for the MCMC sampling
iterations <- 20000
# Number of chains to run
chains <- 4
SCALE <- 1
ns <- iterations*chains/2
# number of sigfigs to show on the table
nfigs <- 4

mod = '1 + TRAITmean + TRAITdiff + STATEmean + STATEdiff'
modelForm = paste('Y ~', mod,'+ (1 | gr(Subj, dist= "student")) + (',mod,'| gr(ROI, dist="student")) + (',mod,'|gr(VOX, dist="student"))')
print('Model Formula:')
print(modelForm)

# get dafualt priors for data and model
priorRBA <- get_prior(formula = modelForm,data=dataset,family = 'student')
# You can assign prior or your choice to any of the parameter in the table below. 
# For example. If you want to assign a student_t(3,0,10) prior to all parameters of class b, 
# the following line does that for you. Parameters in class b are the population effects (cond, STATE and TRAIT)
priorRBA$prior[1] <- "student_t(3, 0, 10)"
priorRBA$prior[6] <- "lkj(2)"
priorRBA$prior[9:11] <- "gamma(3.325,0.1)"
priorRBA$prior[13] <- "gamma(3.325,0.1)"
print("")
# Print the table with priors
print(priorRBA)

# Create a result directory
if (!dir.exists("results")){dir.create("results")}

# Generate the Stan code for our own reference
stan_code <- make_stancode(modelForm,data=dataset,chains = chains,family = 'student',prior = priorRBA)
cat(stan_code,file = "results/stancode.stan",sep='\n')

# Following run the BML model
fm <- brm(modelForm,
          data=dataset,
          chains = chains,
          family = 'student',
          prior = priorRBA,
          inits=0, iter=iterations, 
          control = list(adapt_delta = 0.99, max_treedepth = 15))

# Shows the summary of the model
cat(capture.output(summary(fm)),sep = '\n', append=TRUE)

# Save the results as a RData file
save.image(file="results/results.RData")


