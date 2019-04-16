#---------------
# digoxin.model
#---------------

# Variables
States = {A_central, A_periph}
Inputs = {Dose}
Outputs = {C_central}

# Structural model parameters
k_12 = 1.02;
k_21 = 0.15;
k_10 = 0.18;
V_central = 58.2;

# Measurement error 
Ve_C_central = 1;

# Initalization
Initialize {
  A_central = Dose;
}

# Dynamics
Dynamics {
  # Central compartment quantity
  dt(A_central) = k_21 * A_periph - k_12 * A_central - k_10 * A_central;
  # Peripheral compartment quantity 
  dt(A_periph) = k_12 * A_central - k_21 * A_periph;
}

CalcOutputs { 
	C_central = A_central / V_central ; 
}

End.
