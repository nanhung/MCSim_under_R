# Description ####
# One-compartment model with first-order absorption rate and 
# linear elimination (Based on R httk package)  
#
# version 1 (This version is use to test forward simple pk)
#
# Units: 
# - time in hours
# - volumes in liters
#
# Nan-Hung Hsieh - May 2019

States  = { A_rest };  # Quantity in central compartment (mg)
Outputs = { Ccompartment };  # Concentration in central compartment (uM)
Inputs = { Oral_input }; # Drug input (uM)

# Oral input modeling
Absdose;
IngDose    = 1.0; # ingested input (mg)
Fgutabs    = 1.0; #
Period     = 0.0; # period of the exposure/no exposure cycle (h)
Tlag       = 0.0; # Absorption lagtime (h)
kgutabs    = 0.1; # Intestinal absortion rate (/h)
Oral_input = PerExp (Absdose, Period, Tlag, kgutabs);
# Molecular weight (g/mole)
MW;
# Distribution volumes (L)
Vdist  = 1;
# Elimination rate constants (/h)
kelim = 1;

Initialize{ Absdose = IngDose * Fgutabs; }
Dynamics { dt (A_rest) = kgutabs * Oral_input - kelim * A_rest; }
CalcOutputs { Ccompartment  = A_rest  / Vdist / MW * 1000; }

End.
