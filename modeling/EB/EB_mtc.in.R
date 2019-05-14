Integrate (Lsodes, 1e-9, 1e-11 , 1);

MonteCarlo ("", 10, 95814);

Distrib (BW, Uniform, 0.035, 0.055);

Distrib (Pb, Uniform, 21.35, 85.4);
Distrib (Pl, Uniform, 0.98, 3.92);
Distrib (Pf, Uniform, 18.2, 72.8);
Distrib (Pm, Uniform, 0.3045, 1.218);
Distrib (Pvrg, Uniform, 0.98, 3.92);
Distrib (Ppu, Uniform, 0.98, 3.92);
Distrib (Pbr, Uniform, 0.98, 3.92);

Distrib (VmaxC, Uniform, 1.5975, 25.56);
Distrib (VmaxClu, Uniform, 3.35, 53.6);
Distrib (VmaxCvr, Uniform, 4.35, 69.6);

Simulation { # 3 100 ppm - 4 hr
  # Inhalation concentration in ppm
  Cppm = NDoses (2, 100, 0, 0, 4 ); 
  PrintStep(Cvtot, 0, 6, 0.1);  
} 

END.
