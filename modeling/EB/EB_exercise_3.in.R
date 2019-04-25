# ./mcsim.EB_v2.model.R.exe EB_exercise_3.in.R

Integrate (Lsodes, 1e-9, 1e-11 , 1);

Simulation { # 1 1 ppm
  # Inhalation concentration in ppm
  Cppm = NDoses (2, 1, 0, 0, 8 ); 
  PrintStep(Ain, Amet_Rl, Amet_Rlu, Amet_Rvrg, Amet, 0, 8, 0.5);  
} 

Simulation { # 2 10 ppm
  # Inhalation concentration in ppm
  Cppm = NDoses (2, 10, 0, 0, 8 ); 
  PrintStep(Ain, Amet_Rl, Amet_Rlu, Amet_Rvrg, Amet, 0, 8, 0.5);  
} 

Simulation { # 3 100 ppm
  # Inhalation concentration in ppm
  Cppm = NDoses (2, 100, 0, 0, 8 ); 
  PrintStep(Ain, Amet_Rl, Amet_Rlu, Amet_Rvrg, Amet, 0, 8, 0.5);  
} 

Simulation { # 4 1000 ppm
  # Inhalation concentration in ppm
  Cppm = NDoses (2, 1000, 0, 0, 8 ); 
  PrintStep(Ain, Amet_Rl, Amet_Rlu, Amet_Rvrg, Amet, 0, 8, 0.5);  
} 

END.