# ./mcsim.EB.model.R.exe EB_exercise_1.in.R

Integrate (Lsodes, 1e-9, 1e-11 , 1);

Simulation { 
  # Inhalation concentration in ppm
  Cppm = NDoses (2, 100, 0, 0, 4 ); 
  PrintStep(Cvtot, 0, 6, 0.01);  
} 

END.
