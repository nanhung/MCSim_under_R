#-------------------
# dogoxin.mcmc.in
#-------------------
Integrate (Lsodes, 1e-4, 1e-6, 1);

MCMC ("digoxin.mcmc.out","","",50000,0,10,50000,1111);

Level { # top
  Distrib(k_12, TruncLogNormal, 0.2, 4, 0.001, 2);
  Distrib(k_21, TruncLogNormal, 0.2, 4, 0.001, 2);
  Distrib(k_10, TruncLogNormal, 0.2, 4, 0.001, 2);
  Distrib(V_central, TruncLogNormal, 40, 4, 40, 60);
  Distrib(Ve_C_central_SD, HalfNormal, 4);
  
  Likelihood(C_central , Normal, Prediction(C_central) , Ve_C_central_SD);
  
  Simulation {
  Dose = 509;
  Print (C_central, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 23);
  Data (C_central, 4.6244, 2.7654, 1.3224, 0.9563, 0.8843, 0.8648, 0.8363, 0.7478, 0.7232, 0.5655);
  }
} End. 
