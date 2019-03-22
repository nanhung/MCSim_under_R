## Pre-setting and load functions ####
# clear() # use to clear the exe and out file
# make sure to put the model and input files to "input folder"
source("FUN.R") 


## simple model #####
# Define the name of model and input files
mName <- "simple.model.R" # the model file put in the model folder
inName <- "simple.in.R" # the input file put in the infile folder

# Create the executable file
compile_mod(mName)
# file.exists("mcsim_simple.model.R.exe") # check if you create the "mcsim_simple.model.exe" file successfully

# Run!!
out <- run_mcsim(mName, inName)

# Print result
out

# Plot 
out <- out[2:13,] # omit 0
pkplot(out, var = "y0", log = "x", xlab = "", ylab = "", type = "b")
pkplot(out, var = "y1", log = "x", xlab = "", ylab = "", type = "b")
pkplot(out, var = "y2", log = "x", xlab = "", ylab = "", type = "b")

## Digoxin MCMC ####
# Define the input variable
mName <- "digoxin.model.R" # the model file put in the model folder
inName <- "digoxin.mcmc.in.R" # the input file put in the infile folder

# Create the executable file
compile_mod(mName)

# Run!!
set.seed(1111)
out <- run_mcsim(mName, inName)

par(mfrow= c(2,2))
names(out)
plot(out$k_12.1., type = "l")
plot(out$k_21.1., type = "l")
plot(out$k_10.1., type = "l")
plot(out$V_central.1., type = "l")

check_df <- read.delim("chk.out")
par(mfrow= c(2,1))
plot(check_df$Time, check_df$Data)
lines(check_df$Time, check_df$Prediction)
plot(check_df$Data, check_df$Prediction)
abline(0,1)
