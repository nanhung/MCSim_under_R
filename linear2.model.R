#------------------------------------------------------------------------------
# linear.model
#
# Linear Model with random noise added.
# P = A + B * time + Normal (0, SD_true)
# Setting SD_true to zero gives the deterministic version, which can be used
# as link function in statistical models.
#
# Copyright (c) 1993-2008 Free Software Foundation, Inc.
#------------------------------------------------------------------------------

Outputs = {y,z};

# Model Parameters
A = 0;
B = 1;

# Statistical parameter
SD_y = 1;
SD_z = 1;

CalcOutputs { 
  y = A + B * t; 
  z = A + B * t;
} 

End.
