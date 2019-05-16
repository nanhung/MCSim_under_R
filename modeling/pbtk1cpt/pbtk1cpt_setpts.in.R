# Description ####
## ./mcsim.pbtk1cpt.model.R.exe pbtk1cpt_setpt.in.R 

Integrate (Lsodes, 1e-9, 1e-9, 1);

SetPoints ("", "setpts.out", 0, Vdist, kelim, kgutabs, Fgutabs);

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
