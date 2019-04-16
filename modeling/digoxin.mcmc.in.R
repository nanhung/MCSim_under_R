#-------------------
# dogoxin.mcmc.in
#-------------------
Integrate (Lsodes, 1e-4, 1e-6, 1);

MCMC("sim.out","", # name of output and restart file
     "",           # name of data file
     2000,0,      # iterations, print predictions flag,
     10,2000,     # printing frequency, iters to print
     10101010);    # random seed (default )

Level { # top
  Distrib(k_12, LogUniform, 0.01, 10);
  Distrib(k_21, LogUniform, 0.01, 10);
  Distrib(k_10, LogUniform, 0.01, 10);
  Distrib(V_central, TruncNormal, 50, 5, 40, 60);
  Distrib(Ve_C_central, LogUniform, 0.01, 0.5); # 10% to 70% residual error
  
  Likelihood(C_central , Normal, Prediction(C_central) , Ve_C_central);
  
  Simulation {
  Dose = 509;
  Print (C_central, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 23);
  Data (C_central, 4.6244, 2.7654, 1.3224, 0.9563, 0.8843, 0.8648, 0.8363, 0.7478, 0.7232, 0.5655);
  }
} End. 
