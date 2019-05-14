

SetPoints ("setpoint.out", "sim.out", 0, Ka, Ke);

Simulation { # 1
  OralDose = 100; 
  BW = 60;
  PrintStep (C_central, 0, 8, 1);
}

End.
