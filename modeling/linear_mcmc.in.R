## linear_mcmc.in.R ####
MCMC ("MCMC.default.out","", # name of output and restart file
     "",           # name of data file
     2000,0,       # iterations, print predictions flag,
     1,2000,       # printing frequency, iters to print
     10101010);    # random seed (default )

Level {
  
  Distrib(A, Normal, 0, 2); # prior of intercept coefficient
  Distrib(B, Normal, 1, 2); # prior of slope coefficient
  
  Likelihood(y, Normal, Prediction(y), 0.5); #  # exact SD
  
  Simulation {
    PrintStep (y, 0, 10, 1); #seq(0, 10 1)
    Data  (y, 0.0, 0.15, 2.32, 4.33, 4.61, 6.68, 7.89, 7.13, 7.27, 9.4, 10.0);
  }
}

End.
