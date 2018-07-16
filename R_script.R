# Download all files form this repository (the current version of MCSim is 6.0.1)
Sys.getenv("PATH")

# You have two options to use GNU compiler:
# If you have installed MinGW in your PC you can use
Sys.setenv(PATH = paste("c:\\MinGW\\bin", Sys.getenv("PATH"), sep=";"))

# Otherwise, if you have installed Rtools, you can assign the bin location manually,
# Sys.setenv(PATH = paste("c:\\Rtools\\bin", Sys.getenv("PATH"), sep=";"))
# Sys.setenv(PATH = paste("c:\\Rtools\\mingw_32/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(PATH = paste("c:\\Rtools\\mingw_64/bin", Sys.getenv("PATH"), sep=";"))
# Sys.setenv(BINPREF = "c:\\Rtools\\mingw_64/bin/")

# Check the GNU compiler 
Sys.which("gcc")
system('g++ -v')

# Make sure you are in the correct working directory 
getwd()

# Define the name of model and input files
mName <- "simple.model"
inName = "simple.in"

# Create mod.exe in "mod" folder
system(paste("gcc -o ./mod/mod.exe ./mod/*.c ", sep = "")) 

# Compile the "simple.model" to "simple.c" 
system(paste("./mod/mod.exe ", mName, " ", mName, ".c", sep = "")) 

# Compile the "simple.model.c" to the executable program named "mcsim.simple.model.exe"
system(paste("gcc -O3 -I.. -I./sim -o mcsim_", mName, ".exe ", mName, ".c ./sim/*.c -lm ", sep = ""))

# Run!!
system(paste("./mcsim_", mName, ".exe ", inName, sep = ""))

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
