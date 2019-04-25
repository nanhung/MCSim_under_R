## Description ####
# EBZ.model.R 
#
# Pratice file (Translate from Berkely Madonna to GNU MCSim) 
# A PBPK model for Ethylbenzene toxicokinetics in male rats.
#
# Units:
# - time in h
# - volumes in L
# - massess in moles
# - concentration in moles/L

## States ####
States = {
  
};

##  Outputs ####
Outputs = {
  
};

## Inputs ####
Inputs = { 
  
};  

## Parameters ####
# Exposure Parameter -

# Physiologocal Parameter -

# Substance-specific Parameter -
# molecular weight

# blood/air and tissue/blood partition coefficients

# metabolic parameters

# Scaling Parameter -
# body weight scaling

# conversion factor

# tissue flows, L/hr

# tissue volumes, L

# the Michaelis–Menten kinetics

## Initialize ####
Initialize {

  # body weight scaling 

  # conversion factor

  # tissue flows, L/hr

  # tissue volumes, L

  # the Michaelis–Menten kinetics

}

## Dynamics ####
Dynamics {
  
  # exposure concentration converted to mol/L

  # calculated concentrations of ethylbenzene

  # intake rate

  # metabolism

  # ethylbenzene uptake and metabolism

  # amount of ethylbenzene metabolize

  # amount of ethylbenzene intake

}

## CalcOutputs ####
CalcOutputs {
  
  # concentration

  # mass Balance
  
}

End.
