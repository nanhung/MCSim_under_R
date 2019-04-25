Integrate (Lsodes, 1e-9, 1e-11 , 1);

Simulation { # 1 100 ppm - 8 hr
  # Inhalation concentration in ppm
  Cppm = NDoses (2, 100, 0, 0, 8 ); 
  PrintStep(Cart, Cvipu,	Cvtot, 0, 48, 0.1);  
} 

Simulation { # 2 1 mg/m3 - 8 hr
  # Inhalation concentration in ppm
  Cppm = NDoses (2, 0.00023, 0, 0, 8 ); 
  PrintStep(Cart, Cvipu,	Cvtot, 0, 48, 0.1);  
} 

Simulation { # 3 100 ppm - 4 hr
  # Inhalation concentration in ppm
  Cppm = NDoses (2, 100, 0, 0, 4 ); 
  PrintStep(Cvtot, 0, 6, 0.1);  
} 

END.
