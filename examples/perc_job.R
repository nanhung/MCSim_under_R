message(paste0("Starting time: ", Sys.time()))
source("MCSim/function.R")
system.time(out_2 <- mcsim(model = "perc.model.R", input = "perc_mcmc-hier.in.R", dir = "modeling/perc", parallel = T))
message(paste0("Ending time: ", Sys.time()))