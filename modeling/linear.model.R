## linear.model.R ####
Outputs = {y}

# Model Parameters
A = 0; # 
B = 1;

# Statistical parameter
SD_true = 0;

CalcOutputs { y = A + B * t + NormalRandom(0,SD_true); }

End.
