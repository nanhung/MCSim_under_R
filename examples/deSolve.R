## simple model ####
library(deSolve)
model <- "simple.model.R"
makemcsim(model = model, deSolve = T, dir = "modeling/simple")
times <- c(0, 0.4*10^(0:11)) # define parameter values
Y <- c(y1 = 1.0, y2 = 0.0, y3 = 0.0) # define initial state values
parms <- c(k1 = 0.04, k2 = 1e4, k3=3e7) 
out <- ode(Y, times, func = "derivs", parms = parms, dllname = model, 
           initfunc = "initmod", nout = 1, outnames = "Sum")
plot(out, log="x",type="b")


## digoxin ####
clear()
model <- "digoxin.model.R"
makemcsim(model = model, deSolve = T)
parms <- initParms()
newParms <- c(parms, Dose = 10)
Y <- initStates(newParms)
newParms
Y

times <- seq(0, 24, 1)

out <- ode(Y, times, func = "derivs", parms = parms, dllname = model, 
           initfunc = "initmod", nout = 1, outnames = "C_central")
out

plot(out)


