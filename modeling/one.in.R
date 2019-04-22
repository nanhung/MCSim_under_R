# File name: one.in.R
# ./mcsim.one.modle.R.exe one.in.R

Integrate (Lsodes, 1e-12, 1e-15, 1);

BW      = 70;
Period  = 1E6;  # One-time dose
Ka      = 1.3;  
Pct_M_central = 1;

Simulation {
  PrintStep (Oral_input, A_central, A_elim, A_total, C_central, 0, 96, 0.5)
}

End.
