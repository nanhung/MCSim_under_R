# ./mcsim.perc.model.R.exe perc.lsodes.in.R
Integrate (Lsodes, 1e-6, 1e-8, 1);

Simulation {
  
  InhMag = 72;            # ppm
  Period = 1e10;          # Only one dose
  Exposure = 240;         # 4 hour exposure
  
  # measurements before end of exposure
  # and at [5' 30'] 2hr 18 42 67 91 139 163
  
  Print (C_exh_ug, C_ven, Pct_metabolized, 239.9 245 270 360 1320 2760 4260 5700 8580 10020 );
  
}

Simulation {
  
  InhMag = 144;            # ppm
  Period = 1e10;          # Only one dose
  Exposure = 240;         # 4 hour exposure
  
  # measurements before end of exposure
  # and at [5' 30'] 2hr 18 42 67 91 139 163
  
  Print (C_exh_ug, C_ven, Pct_metabolized, 239.9 245 270 360 1320 2760 4260 5700 8580 10020 );
  
}

END.
