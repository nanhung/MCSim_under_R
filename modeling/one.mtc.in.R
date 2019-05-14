## ./mcsim.one.model.R.exe one.mtc.in.R ####
## Monte Carlo simulation input file for one compartment model

MonteCarlo ("", 10, 95814);

Distrib (Ka,  Uniform, 0.2,    0.8);
Distrib (Ke,  Uniform, 0.03,   0.1);

Simulation { # 1
  OralDose = 100; 
  BW = 60;
  PrintStep (C_central, 0, 8, 1);
}

End.
