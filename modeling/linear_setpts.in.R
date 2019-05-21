## ./mcsim.linear.model.R.exe linear_setpt.in.R ####
SetPoints ("", "setpts.out", 0, A, B);

Simulation { 
  PrintStep (y, 0, 10, 1);
  Data  (y, -0.0289654, 1.15968, 2.32502, 3.33289, 4.61105, 5.6818, 
         6.89044, 8.13242, 9.27033, 10.4522, 11.6703);
} 

END.
