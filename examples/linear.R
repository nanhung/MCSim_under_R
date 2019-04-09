## Pre-setting and load functions ####
source("function.R") 

## Linear modeling ####
mName <- "linear.model.R" # model file
inName <- "linear.in.R" # input file; simple simulation w/ given parameters

# Create the executable file 'mcsim.linear.model.R.exe'
makemcsim(mName)

# Run!!
out <- mcsim(mName, inName) #'./mcsim.linear.model.R.exe linear.in.R'
out # check result

i = 1:31
# i = 34:36
x <- as.numeric(as.character(out$Time[i]))
y <- as.numeric(as.character(out$y[i]))
plot(x, y, type = "b") # 

## Linear mcmc ####
# Assign the input file
inName <- "linear.mcmc.in.R"

# Generate the first chain
set.seed(1111) 
out <- mcsim(mName, inName) #'./mcsim.linear.model.R.exe linear.mcmc.in.R'
out # check result

plot(out$A.1., type = "l")
plot(out$B.1., type = "l")

# Density plot
plot(density(out$A.1.[i]))
plot(density(out$B.1.[i]))

# Visualizing the fitting result 
chk <- read.delim("chk.out")
chk
plot(chk$Time, chk$Data)
lines(chk$Time, chk$Prediction)

# Check convergence
set.seed(2234) 
out2 <- mcsim(mName, inName) # Generate the 2nd chain
set.seed(3234) 
out3 <- mcsim(mName, inName) # Generate the 3rd chain
set.seed(4234) 
out4 <- mcsim(mName, inName) # Generate the 4th chain

# Use rstan package to check result
library(rstan)
n.iters <- nrow(out)
n.chains <- 4
n.parms <- 2
parms_name <- c("A.1.","B.1.")
sims = array(0, c(n.iters, n.chains, n.parms))
sims[,1,] = as.matrix(out[, parms_name])
sims[,2,] = as.matrix(out2[, parms_name])
sims[,3,] = as.matrix(out3[, parms_name])
sims[,4,] = as.matrix(out4[, parms_name])
report <- monitor(sims, digit=4)
row.names(report) <- parms_name
report


