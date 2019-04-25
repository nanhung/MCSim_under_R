## Description ####
# EB.model.R
#
# A PBPK model for Ethylbenzene toxicokinetics in male rats.
# 
# Version: 1
# 
# Units:
# - time in h
# - volumes in L
# - massess in moles
# - concentration in moles/L

## States ####
States = {
  Af,         # Amount of ethylbenzene in fat (moles)
  Al,         # Amount of ethylbenzene in liver
  Am,         # Amount of ethylbenzene in muscle
  Avrg,       # Amount of ethylbenzene in richly perfused tissues
  Apu,        # Amount of pulmonary ethylbenzene
  Abr,        # Amount of bronchial ethylbenzene
  AUCvtot,    # AUC of blood ethylbenzene (mol-hr/L)
  Ain,        # Amount of inhaled ethylbenzene
  Amet        # Amount metabolized (moles)
};

##  Outputs ####
Outputs = {
  Cart,    # Concentration in arterial blood (mol/L)
  Cvipu,   # ~             in pulmonary
  Cvtot,   # ~             in venous blood
  At       # ~             in all tissues
};

## Inputs ####
Inputs = { Cppm };    # Concentration (inhaled)

## Parameters ####
# Exposure parameter
Cppm = 75;

# Physiologocal parameter
BW = 0.043;  # male rat BW

# substance-specific parameter
# molecular weight
MW = 106.16;

# blood/air and tissue/blood partition coefficients
Pb = 42.7;
Pl = 1.96;
Pf = 36.4;
Pm = 0.609;
Pvrg = 1.96;
Ppu = 1.96;
Pbr = 1.96;

# metabolic parameters
VmaxC = 6.39;     #liver
VmaxClu = 13.4;   #lung
VmaxCvr = 17.4;   #richly perfused

# Scaling parameter - 
# body weight scaling
BW74;

# conversion factor
Cfac;   # mg to mol

# tissue flows, L/hr
Qtot;   # total
Qalv;   # alveolar
Qpu;    # pulmonary
Qbr;    # bronchial
Qf;     # fat
Qvrg;   # richly perfused
Ql;     # liver
Qm;     # muscle

# tissue volumes, L
Vf;     # fat
Vl;     # liver
Vm;     # muscle
Vvrg;   # richly perfused
Vlu;    # lung
Vpu;    # pulmonary
Vbr;    # bronchial

# the Michaelis–Menten kinetics
Vmax;
Km;
Vmaxlu;
Kmlu;
Vmaxvr;
Kmvr;

## Initialize ####
Initialize {
  # Body weight scaling
  BW74 = pow(BW, 0.74);
  
  # conversion factor
  Cfac = MW * 1000;    # mg to mol
  
  # tissue flows, L/hr
  Qtot = 15 * BW74;
  Qalv = 15 * BW74;
  Qpu = 0.928 * Qtot;
  Qbr = 0.072 * Qtot;
  Qf = 0.09 * Qpu;
  Qvrg = 0.51 * Qpu;
  Ql = 0.25 * Qpu;
  Qm = 0.15 * Qpu;
  
  # tissue volumes, L
  Vf = 0.06 * BW;       #fat
  Vl = 0.04 * BW;       #liver
  Vm = 0.76 * BW;       #muscle
  Vvrg = 0.05 * BW;     #richly perfused
  Vlu = 0.014 * BW;     #lung
  Vpu = 0.454 * Vlu;    #pulmonary
  Vbr = 0.545 * Vlu;    #bronchial
  
  # the Michaelis–Menten kinetics
  Vmax = VmaxC * BW74 / Cfac;
  Km = 1.04 / Cfac;
  Vmaxlu = VmaxClu * BW74 / Cfac;
  Kmlu = 5.57 / Cfac;
  Vmaxvr = VmaxCvr * BW74 / Cfac;
  Kmvr = 2.33 / Cfac;
}

## Dynamics ####
Dynamics {
  # exposure concentration converted to mol/L
  Cmpl = Cppm * 1E-6 / 24.45;
  Cair = Cmpl;
  
  # calculated concentrations of ethylbenzene
  Cvpu = Apu / (Vpu * Ppu);
  Cvbr = Abr / (Vbr * Pbr);
  Cart = (Qpu * Cvpu + Qbr * Cvbr) / Qtot;
  Cvf = Af / (Vf * Pf);
  Cvl = Al / (Vl * Pl);
  Cvvrg = Avrg / (Vvrg * Pvrg);
  Cvm = Am / (Vm * Pm);
  Cvtot = (Ql * Cvl + Qf * Cvf + Qm * Cvm + Qvrg * Cvvrg) / Qpu;
  Cvipu = Pb * (Qalv * Cair + Qpu * Cvtot) / (Pb * Qpu + Qalv);
  
  # intake rate
  Calv = Cvipu / Pb;
  RI = Qalv * (Cair-Calv);
  
  # metabolism
  Rlu = Vmaxlu * Cvbr / (Kmlu + Cvbr);
  Rl = Vmax*Cvl / (Km + Cvl);
  Rvrg = Vmaxvr * Cvvrg / (Kmvr + Cvvrg);
  
  # ethylbenzene uptake and metabolism
  dt(Apu) = Qpu * (Cvipu - Cvpu);
  dt(Abr) = Qbr * (Cart - Cvbr) - Rlu;
  dt(Al) = Ql * (Cart - Cvl) - Rl;
  dt(Af) = Qf * (Cart - Cvf);
  dt(Avrg) = Qvrg * (Cart - Cvvrg) - Rvrg;
  dt(Am) = Qm * (Cart - Cvm);
  dt(AUCvtot) = Cvtot;
  
  # amount of ethylbenzene metabolize
  dt(Amet) = Rl + Rlu + Rvrg;
  
  # mass Balance
  dt(Ain) = RI;
  At = Abr + Apu + Al + Af + Avrg + Am;
}

## CalcOutputs ####
CalcOutputs {
  
  # Concentration
  Cart = (Cart < 0 ? 0 : Cart);
  Cvtot = (Cvtot < 0 ? 0 : Cvtot);
  Cvipu = (Cvipu < 0 ? 0 : Cvipu);
  
  # Mass Balance
  MBal = (1 + At + Amet) / (1 + Ain);
}

End.
