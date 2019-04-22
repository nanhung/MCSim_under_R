## Description ####
# 1-compartment model with 1st-order absorption rate and 
# linear elimination 
#
# version: 1
#
# Date: 04-25-2019
# 
# Units: 
# - time in hr
# - volumes in L
# - masses of substances in mg
# - concentrations of substances in mg/L

## States ####
States = { A_central,  # Quantity in central compartment (mg)
           A_elim}     # ~        eliminated

## Outputs ####
Outputs = {C_central,  # Concentration in central compartment (mg/L)
           A_total}    # Total quantity for mass balance

## Inputs ####
Inputs = {Oral_input} # Chemical input (mg)

## Parameters ####

# Chemical-specific parameters
Ke = 0.1;             # Elimination rate constant (1/h)
Pct_M_central = 0.05; # % body weight

# Physiological-specific parameter
BW = 60; # Body weight (kg) 

# Exposre parameter
OralDose   = 100;  # Oral dose (mg/kg)
Period     = 12.0; # period of the exposure/no exposure cycle (h)
Tlag       = 0.0;  # Absorption lagtime (h)
Ka         = 0.1;  # Intestinal absorption rate constant (1/h)

# Scale parameter computed in Initialize
V_central; # Distribution volume of central compartment (L)
IngDose;   # Ingested dose (mg)

Oral_input = PerExp (IngDose, Period, Tlag, Ka);

## Initialize ####
Initialize {  
  IngDose = BW * OralDose;
  V_central = BW * Pct_M_central; 
}

## Dynamics ####
Dynamics {
  dt (A_elim)    = Ke * A_central;
  dt (A_central) = Ka * Oral_input - dt(A_elim);
}

## CalcOutputs ####
CalcOutputs {
  C_central = A_central / V_central;
  A_total   = A_central + A_elim;
}

End.
