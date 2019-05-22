## ./mcsim.linear.model.R.exe linear_setpt.in.R ####
SetPoints ("", "setpts.out", 0, A, B);

Simulation { 
  PrintStep (y, 0, 10, 1);
  Data  (y, 0.0, 0.15, 2.32, 4.33, 4.61, 6.68, 7.89, 7.13, 7.27, 9.4, 10.0);
} 

END.
