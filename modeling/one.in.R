# ./mcsim.one.model.R.exe one.in.R

Integrate (Lsodes, 1e-12, 1e-15, 1);

Period  = 1E6;  # One-time dose
Ka      = 1.3;  
Pct_M_central = 1;

Simulation { # 1
  OralDose = 100; 
  BW = 60;
  PrintStep (Oral_input, A_central, A_elim, A_total, C_central, 0, 24, 1);
}

Simulation { # 2
  OralDose = 150;
  BW = 80;
  PrintStep (Oral_input, A_central, A_elim, A_total, C_central, 0, 24, 1);
}

End.
