# Description ####
## ./mcsim.pbtk1cpt.model.R.exe pbtk1cpt.mtc.in.R 

Integrate (Lsodes, 1e-06, 1e-06 , 1);

MonteCarlo ("",1000,11920);

Distrib ( Vdist,Normal_cv,6.137241,0.2);
Distrib ( kelim,Normal_cv,0.02283233,0.2);
Distrib ( kgutabs,Normal_cv,2.18,0.2);
Distrib ( Fgutabs,Uniform,0.8,1);

#---------------------------------------- 
# Simulation scenario
#----------------------------------------

Simulation { 
  
  MW = 228.291;
  Period = 24;
  IngDose = 1.0;
  
  PrintStep (Ccompartment, 0, 23, 1);
} 
END.
