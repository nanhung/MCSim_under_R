## ./mcsim.one.model.R.exe one_setpts.in.R ####
## Setpoint simulation input file for one compartment model
## Use "sim.out" that generated from "one_mtc.in.R"

SetPoints ("setpts.out", "simmc.out", 0, Ka, Ke);

Simulation { # 1
  OralDose = 100; 
  BW = 60;
  PrintStep (C_central, 0, 8, 1);
}

End.
