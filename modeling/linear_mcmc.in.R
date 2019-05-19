#-------------------
# linear_mcmc.in
#-------------------

MCMC ("MCMC.default.out","", # name of output and restart file
     "",           # name of data file
     4000,0,       # iterations, print predictions flag,
     1,4000,       # printing frequency, iters to print
     10101010);    # random seed (default )

Level {
  
  Distrib(A, Normal, 0, 2); # prior of intercept coefficient
  Distrib(B, Normal, 1, 2); # prior of slope coefficient
  
  Likelihood(y, Normal, Prediction(y), 0.05); #  # exact SD
  
  Simulation {
    PrintStep (y, 0, 10, 1); #seq(0, 10 1)
    Data  (y, -0.0289654, 1.15968, 2.32502, 3.33289, 4.61105, 5.6818, 
           6.89044, 8.13242, 9.27033, 10.4522, 11.6703);
  }
}


End.
  