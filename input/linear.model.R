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

Outputs = {y};

# Model Parameters
A = 0;
B = 1;
SD_true = 0;

# Statistical parameter
Sigma2 = 1;

CalcOutputs { y = A + B * t + NormalRandom(0,SD_true); } 

End.
