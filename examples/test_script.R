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







