#------------------------------------------------------------------------------
# linear.in
#
# Linear model, simple input file
#
# Copyright (c) 1993-2008 Free Software Foundation, Inc.
#------------------------------------------------------------------------------

Simulation {
  
  A  = 0;
  B  = 1;
  SD_true = 0.05;
  
  PrintStep (y, 0, 30, 1);
}

Simulation {
  
  A  = 1;
  B  = 2;
  SD_true = 0.05;
  
  Print (y, 0, 1, 2);
}

END.
