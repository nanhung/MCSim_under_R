message(paste0("Starting time: ", Sys.time()))
source("MCSim/function.R")
system.time(out <- mcsim(model = "EB.model.R", input = "EB_mcmc.in.R", dir = "modeling/EB", parallel = T))
message(paste0("Ending time: ", Sys.time()))