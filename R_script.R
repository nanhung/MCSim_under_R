## Pre-setting and load functions ####
# clear() # use to clear the exe and out file
# make sure to put the model and input files to "input folder"
source("function.R") 

## Linear regression ####
mName <- "linear.model.R" 
inName <- "linear.in.R" # simple simulation w/ given parameters

# Create the executable file 'mcsim.linear.model.R.exe'
makemcsim(mName)

# Run!!
out <- mcsim(mName, inName) #'./mcsim.linear.model.R.exe linear.in.R'
out

set.seed(1111)
inName <- "linear.mcmc.in.R" 


## Digoxin MCMC ####
# Define the input variable
mName <- "digoxin.model.R" # the model file put in the model folder
inName <- "digoxin.mcmc.in.R" # the input file put in the infile folder

# Create the executable file
makemcsim(mName)

# Run!!
set.seed(1111)
out <- mcsim(mName, inName)
out

par(mfrow= c(2,2))
names(out)
plot(out$k_12.1., type = "l")
plot(out$k_21.1., type = "l")
plot(out$k_10.1., type = "l")
plot(out$V_central.1., type = "l")

check_df <- read.delim("chk.out")
par(mfrow= c(1,1))
plot(check_df$Time, check_df$Data)
lines(check_df$Time, check_df$Prediction)
plot(check_df$Data, check_df$Prediction)
abline(0,1)


## simple model #####
# Define the name of model and input files
mName <- "simple.model.R" # the model file put in the model folder
inName <- "simple.in.R" # the input file put in the infile folder

# Create the executable file
makemcsim(mName)
# file.exists("mcsim_simple.model.R.exe") # check if you create the "mcsim_simple.model.exe" file successfully

# Run!!
out <- mcsim(mName, inName)

# Print result
out

# Plot 
out <- out[2:13,] # omit 0
par(mfrow=c(2,2))
plot(out$Time, out$y0, log = "x", type = "b")
plot(out$Time, out$y1, log = "x", type = "b")
plot(out$Time, out$y2, log = "x", type = "b")







