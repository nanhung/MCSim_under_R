# Download all files form this repository (the current version of MCSim is 5.6.5)
Sys.getenv("PATH")

# You have two options to use GNU compiler:
# If you have installed MinGW you can use (I use MinGW in this case)
Sys.setenv(PATH = paste("c:\\MinGW\\bin", Sys.getenv("PATH"), sep=";"))

# If you have installed Rtools you can use
Sys.setenv(PATH = paste("c:\\Rtools\\bin", Sys.getenv("PATH"), sep=";"))


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
par(mfrow = c(2,2))
plot(out[,1], out.3[,2], log="x",type="b", xlab = "time", ylab = "", main = "y1")
plot(out[,1], out.3[,3], log="x",type="b", xlab = "time", ylab = "", main = "y2")
plot(out[,1], out.3[,4], log="x",type="b", xlab = "time", ylab = "", main = "y3")
