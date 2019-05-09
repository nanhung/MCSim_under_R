Outputs = {Dose}

# Model Parameters
Fgutabs = 1; 
kelim = 1;
Vdist = 1;
Css = 1;

CalcOutputs { Dose = Css * kelim * Vdist / Fgutabs * t  ; }

End.