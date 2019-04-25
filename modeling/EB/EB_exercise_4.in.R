# ./mcsim.EB_v2.model.R.exe EB_exercise_4.in.R

Integrate (Lsodes, 1e-9, 1e-11 , 1);

Simulation { # 1 
  Dmgkg = PerExp(180, 1e2, 0.0, 1.0);
  PrintStep(Dmgkg, Cvtot, 0, 24, 0.1);  
} 

End.
