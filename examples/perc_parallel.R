library(foreach)     
library(doParallel)  

detectCores()
cores <- 3
cl <- makeCluster(cores)
registerDoParallel(cl)

# Create executable program
makemcsim("perc.model.R", dir = "modeling/perc")

# Parallel computing
strt<-Sys.time()
system.time( 
  out <- foreach(i = 1:3) %dopar% { mcsim(model = "perc.model.R", input = "perc_mcmc.in.R", dir = "modeling/perc", parallel = T)  }
)
print(Sys.time()-strt)

stopCluster(cl)   
