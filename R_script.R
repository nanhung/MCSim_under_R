# Load functions
source("functions.R")

# Download all files form this repository (the current version of MCSim is 6.0.1)
# Check the chek whether the compiler is in the PATH by using
Sys.getenv("PATH") 

# You have two options to use GNU compiler:
# If you have installed MinGW in your PC you can use
Sys.setenv(PATH = paste("c:\\MinGW\\bin", Sys.getenv("PATH"), sep=";"))

# Otherwise, if you have Rtools installed, you can assign the bin location manually,
# Sys.setenv(PATH = paste("c:\\Rtools\\bin", Sys.getenv("PATH"), sep=";"))
# Sys.setenv(PATH = paste("c:\\Rtools\\mingw_32/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(PATH = paste("c:\\Rtools\\mingw_64/bin", Sys.getenv("PATH"), sep=";"))
# Sys.setenv(BINPREF = "c:\\Rtools\\mingw_64/bin/") # Danger zone

# Check the GNU compiler 
Sys.which("gcc")
system('g++ -v')

# Make sure you are in the correct working directory 
getwd()

# Create mod.exe in "mod" folder
compile_mcsim()
# file.exists("mod/mod.exe") # check if you compile the "mod.exe" file successfully


## simple model #####
# Define the name of model and input files
mName <- "simple.model.R" # the model file put in the model folder
inName <- "simple.in.R" # the input file put in the infile folder

# Create the executable file
compile_mod(mName)
# file.exists("mcsim_simple.model.R.exe") # check if you create the "mcsim_simple.model.exe" file successfully

# Run!!
run_mcsim(mName, inName)

# load output
out <- read.delim("sim.out", skip = 2)

# Print
head (out)

# Plot 
out <- out[2:13,] # omit 0
par(mfrow = c(2,2))
plot(out[,1], out[,2], type="b", log = "x", xlab = "time", ylab = "", main = "y1")
plot(out[,1], out[,3], type="b", log = "x", xlab = "time", ylab = "", main = "y2")
plot(out[,1], out[,4], type="b", log = "x", xlab = "time", ylab = "", main = "y3")

## Digoxin MCMC ####
# Define the input variable
mName <- "digoxin.model.R" # the model file put in the model folder
inName <- "digoxin.mcmc.in.R" # the input file put in the infile folder

# Create the executable file
compile_mod(mName)

# Run!!
run_mcsim(mName, inName)

out <- read.delim("digoxin.mcmc1.out")
par(mfrow= c(2,2))
names(out)
plot(out$k_12.1., type = "l")
plot(out$k_21.1., type = "l")
plot(out$k_10.1., type = "l")
plot(out$V_central.1., type = "l")


