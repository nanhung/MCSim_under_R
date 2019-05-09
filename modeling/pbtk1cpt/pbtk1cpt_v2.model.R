## Description ####
# pbtk1cpt.model.R 

# States and Outputs
States = { Aelimination, Agutlument, Ametabolized, Acompartment, AUC};
Outputs = { Ccompartment };

# Parameters
vdist = 0.5;
ke = 0.2;
km = 0.5;
kgutabs = 2.0;

Dynamics {
  Ccompartment = Acompartment / vdist;
  dt (Aelimination) = ke * Agutlument;
  dt (Agutlument) = - kgutabs * Agutlument - dt (Aelimination);
  dt (Ametabolized) = km * Acompartment;
  dt (Acompartment) = kgutabs * Agutlument - dt (Ametabolized);
  dt (AUC) = Ccompartment;
}
End.