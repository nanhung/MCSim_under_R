# Create mod.exe in MCSim (only need to do when the mod.exe doesn;t exist)
# makemod() 

library(rstan)
library(bayesplot)

## Linear mcmc ####
# Assign the input file
model <- "linear.model.R"
inName <- "linear_mcmc.in.R"

# Generate the first chain
set.seed(1111) 
out <- mcsim(model, inName) #'./mcsim.linear.model.R.exe linear.mcmc.in.R'
out # check result

plot(out$A.1., type = "l")
plot(out$B.1., type = "l")

plot(out$A.1., out$B.1., type = "b")

# Density plot
i <- c(ceiling(nrow(out)/2):nrow(out))
plot(density(out$A.1.[i]))
plot(density(out$B.1.[i]))

plot(out$A.1.[i], out$B.1.[i], type = "b")
cor(out$A.1.[i], out$B.1.[i])


# Visualizing the fitting result 
chk <- read.delim("MCMC.check.out")
chk
plot(chk$Time, chk$Data)
lines(chk$Time, chk$Prediction)

# Check convergence
set.seed(2234) 
out2 <- mcsim(model, inName) # Generate the 2nd chain
set.seed(3234) 
out3 <- mcsim(model, inName) # Generate the 3rd chain
set.seed(4234) 
out4 <- mcsim(model, inName) # Generate the 4th chain

sims <- mcmc_array(data = list(out,out2,out3,out4))

parms_name <- c("A.1.","B.1.")


mcmc_trace(sims, pars = parms_name, facet_args = list(ncol = 1, strip.position = "left"))
mcmc_dens_overlay(x = sims[i,,], pars = parms_name)
mcmc_dens_overlay(x = sims[i,,], pars = "LnData")

mcmc_pairs(sims[i,,], pars = parms_name, off_diag_fun = "hex")

monitor(sims[,,parms_name], digit=4) 



