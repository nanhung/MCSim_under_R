#-------------------------------------------------------------------------------
# Simulation file for one-compartment model with first order absorption.
# 
# Units: 
# - time in hours
# - volumes in liters
# - masses of substances in micromoles
# - concentrations of substances in microM
#
#-------------------------------------------------------------------------------

MW = 228.291; # 

Fgutabs    = 0.860811;
IngDose    = 1.0;   # ingested dose (mg)
Period     = 24;     # period of the exposure/no exposure cycle (h)
kgutabs    = 2.18;  # absortion rate (/h)

# Volumes (L)
Vdist = 6.137241; 

# Rate constant (/h)
kelim = 0.02283233;

Simulation {
  PrintStep (Ccompartment, 0, 480, 1); # h
}
End.
