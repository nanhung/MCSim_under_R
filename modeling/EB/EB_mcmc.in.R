#-------------------
# EB_mcmc.in.R
#-------------------
Integrate (Lsodes, 1e-9, 1e-11 , 1);

MCMC ("MCMC.default.out","", # name of output and restart file
     "",           # name of data file
     10000,0,      # iterations, print predictions flag,
     1,10000,     # printing frequency, iters to print
     10101010);    # random seed (default )

Level { # top

  Distrib ( BW, TruncNormal, 0.043, 0.004, 0.035, 0.051);
  Distrib ( Pb, TruncLogNormal, 42.7, 1.6, 21.35, 85.4);
  Distrib ( Pl, TruncLogNormal, 1.96, 1.6, 0.98, 3.92);
  Distrib ( Pf, TruncLogNormal, 36.4, 10, 0.1, 100);
  Distrib ( Pm, TruncLogNormal, 0.609, 1.6, 0.3045, 1.218);
  Distrib ( Pvrg, TruncLogNormal, 1.96, 1.6, 0.98, 3.92);
  Distrib ( Ppu, TruncLogNormal, 1.96, 1.6, 0.98, 3.92);
  Distrib ( Pbr, TruncLogNormal, 1.96, 1.6, 0.98, 3.92);
  Distrib ( VmaxC, TruncLogNormal, 6.39, 10, 0.0639, 639);
  Distrib ( VmaxClu, TruncLogNormal, 13.4, 1.6, 3.35, 53.6);
  Distrib ( VmaxCvr, TruncLogNormal, 17.4, 10, 0.174, 1740);
  
  Likelihood(Cvtot, LogNormal, Prediction(Cvtot) , 1.1);
  
  Simulation {
    Cppm = NDoses (2, 100, 0, 0, 4 );
    
    Print (Cvtot, 4, 4.5, 5, 5.5, 6);
    Data (Cvtot, 1.818011e-05, 1.215147e-05, 8.195177e-06, 5.180859e-06, 3.579503e-06);
  }
} End. 
