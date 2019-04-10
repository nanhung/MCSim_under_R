#------------------------------------------------------------------------------
# linear.in
#
# Linear model, simple input file
#
# Copyright (c) 1993-2008 Free Software Foundation, Inc.
#------------------------------------------------------------------------------

Simulation {
  
  A  = 1; # given value of intercept 
  B  = 2; # given value of slope 
  SD_true = 2; # given SD of noise 
  
  PrintStep (y, 0, 10, 1); 
}

END.
