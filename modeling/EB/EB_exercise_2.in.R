# ./mcsim.EB.model.R.exe EB_exercise_2.in.R

Integrate (Lsodes, 1e-9, 1e-11 , 1);

Simulation { # 1 1 ppm
  # Inhalation concentration in ppm
  Cppm = NDoses (2, 1, 0, 0, 96 ); 
  PrintStep(Cart, Cvtot, 0, 96, 1);  
} 

Simulation { # 2 10 ppm
  # Inhalation concentration in ppm
  Cppm = NDoses (2, 10, 0, 0, 96 ); 
  PrintStep(Cart, Cvtot, 0, 96, 1);  
} 

Simulation { # 3 100 ppm
  # Inhalation concentration in ppm
  Cppm = NDoses (2, 100, 0, 0, 96 ); 
  PrintStep(Cart, Cvtot, 0, 96, 1);  
} 

Simulation { # 4 1000 ppm
  # Inhalation concentration in ppm
  Cppm = NDoses (2, 1000, 0, 0, 96 ); 
  PrintStep(Cart, Cvtot, 0, 96, 1);  
} 

END.