#------------------------------------------------------------------------------
# Rat.EB.model
#
# A PBPK model for Ethylbenzene toxicokinetics in male rats.
#------------------------------------------------------------------------------
#******************************************************************************
#***                  State Variable Specifications                         ***
#******************************************************************************
States = {	
  # ethylbenzene (moles)
  Af, 		# Amount of ethylbenzene in fat
  Al,			# AUC of ethylbenzene in liver
  Am,			# AUC of ethylbenzene in muscle
  Avrg,		# Amount of ethylbenzene in richly perfused
  Apu,		# Amount of pulmonary ethylbenzene
  Abr,		# Amount of bronchial ethylbenzene
  AUCvtot,	# AUC of ethylbenzene in blood (mol-hr/L)
  Ain,		# Amount of inhaled ethylbenzene
  # moles, metabolized
  Amet # Amount metabolized (moles)
};	
#******************************************************************************
#***                  Input Variable Specifications                         ***
#******************************************************************************
Inputs = { 	Cppm	};	# exposure in ppm 
#******************************************************************************
#***                  Output Variable Specifications                        ***
#******************************************************************************
Outputs = {
  Cart,	# Concentration of ethylbenzene in arterial blood (mol/L) 
  # Cf,	
  # Cl,
  # Cvrg,
  # Cm,
  # Cpu,
  # Cbr,
  Cvipu,
  Cvtot,
  At
  # MBal
};		
#******************************************************************************
#***                  Defaults for input parameters                         ***
#******************************************************************************
##-- exposure, ppm
Cppm = 0;
# Model parameters
# -----------------
BW	= 0.043;  # Species-specific defaults during initialization #male rat BW
BW74 = pow(BW, 0.74);
#tissue flows, L/hr
Qtot;
Qalv;
Qpu;	#pulmonary
Qbr;	#bronchial
Qf;		#fat
Qvrg;	#richly perfused
Ql;		#live
Qm;		#muscle
#tissue volumes, L
Vf;		#fat
Vl;		#liver
Vm;		#muscle
Vvrg;		#richly perfused
Vlu;		#lung
Vpu;	#pulmonary
Vbr;	#bronchial
MW = 106.16;
Cfac;	#mg to mol conversion factor
#blood/air and tissue/blood partition coefficients
Pb = 42.7;
Pl = 1.96;
Pf = 36.4;
Pm = 0.609;
Pvrg = 1.96;
Ppu = 1.96;
Pbr = 1.96;
Ain = 0;
#ethylbenzene metabolic parameters, CLh, Vmax mol/hr, Km, M
VmaxC = 6.39;		#liver
VmaxClu = 13.4;		#lung
VmaxCvr = 17.4;		#richly perfused
Vmax;
Km;
Vmaxlu;
Kmlu;
Vmaxvr;
Kmvr;

#******************************************************************************
#***                  Parameter Initialization and Scaling                  ***
#******************************************************************************
Initialize {
  #tissue flows, L/hr
  Qtot = 15*BW74;
  Qalv = 15*BW74;
  Qpu = 0.928*Qtot;
  Qbr = 0.072*Qtot;
  Qf = 0.09*Qpu;
  Qvrg = 0.51*Qpu;
  Ql = 0.25*Qpu;
  Qm = 0.15*Qpu;
  #tissue volumes, L
  Vf = 0.06*BW;		#fat
  Vl = 0.04*BW;		#liver
  Vm = 0.76*BW;		#muscle
  Vvrg = 0.05*BW;		#richly perfused
  Vlu = 0.014*BW;		#lung
  Vpu = 0.454*Vlu;	#pulmonary
  Vbr = 0.545*Vlu;	#bronchial
  Cfac = MW*1000;	#mg to mol conversion factor
  #ethylbenzene metabolic parameters
  Vmax = VmaxC*BW74/Cfac;
  Km = 1.04/Cfac;
  Vmaxlu = VmaxClu*BW74/Cfac;
  Kmlu = 5.57/Cfac;
  Vmaxvr = VmaxCvr*BW74/Cfac;
  Kmvr = 2.33/Cfac;
};

Dynamics {
  # Exposure concentration converted to mol/L
  Cmpl = Cppm*1E-6/24.45;	
  Cair = Cmpl; 
  #calculated concentrations of ethylbenzene
  Cvpu = Apu/(Vpu*Ppu);
  Cvbr = Abr/(Vbr*Pbr);
  Cart = (Qpu*Cvpu + Qbr*Cvbr)/Qtot;
  Cvf = Af/(Vf*Pf);
  Cvl = Al/(Vl*Pl);
  Cvvrg = Avrg/(Vvrg*Pvrg);
  Cvm = Am/(Vm*Pm);
  Cvtot = (Ql*Cvl + Qf*Cvf + Qm*Cvm + Qvrg*Cvvrg)/Qpu;
  Cvipu = Pb*(Qalv*Cair + Qpu*Cvtot)/(Pb*Qpu + Qalv);
  #Intake rate
  Calv = Cvipu/Pb;
  RI = Qalv*(Cair-Calv);
  #metabolism
  Rlu = Vmaxlu*Cvbr/(Kmlu + Cvbr);
  Rl = Vmax*Cvl/(Km + Cvl);
  Rvrg = Vmaxvr*Cvvrg/(Kmvr + Cvvrg);
  #ethylbenzene uptake and metabolism	
  dt(Apu) = Qpu*(Cvipu - Cvpu);
  dt(Abr) = Qbr*(Cart - Cvbr) - Rlu;
  dt(Al) = Ql*(Cart - Cvl) - Rl;
  dt(Af) = Qf*(Cart - Cvf);
  dt(Avrg) = Qvrg*(Cart - Cvvrg) - Rvrg;
  dt(Am) = Qm*(Cart - Cvm);
  dt(AUCvtot) = Cvtot;
  #amount of ethylbenzene metabolize
  dt(Amet) = Rl + Rlu + Rvrg;
  #Mass Balance
  dt(Ain) = RI;
  At = Abr+Apu+Al+Af+Avrg+Am;
};

CalcOutputs {
  Cart = (Cart < 0 ? 0 : Cart);
  # Cf = (Af/Vf < 0 ? 0 : Af/Vf);
  # Cl = (Al/Vl < 0 ? 0 : Al/Vl);
  # Cvrg = (Avrg/Vvrg < 0 ? 0 : Avrg/Vvrg);
  # Cm = (Am/Vm < 0 ? 0 : Am/Vm);
  # Cpu = (Apu/Vpu < 0 ? 0 : Apu/Vpu);
  # Cbr = (Abr/Vbr < 0 ? 0 : Abr/Vbr);
  Cvtot = (Cvtot < 0 ? 0 : Cvtot);
  Cvipu = (Cvipu < 0 ? 0 : Cvipu);
  #Mass Balance
  MBal = (1+At+Amet)/(1+Ain);
};
End.
