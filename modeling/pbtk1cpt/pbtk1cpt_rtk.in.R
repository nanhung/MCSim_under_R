MonteCarlo ("", 1000, 10101010); # Use RandomSeed = 10101010 to apply set.seed() in R

Distrib (Css, LogNormal, 0.00018, 7);
Distrib (Fgutabs, Uniform, 0.8, 1.0);
Distrib (kelim, Normal_cv, 0.02283233, 0.2);
Distrib (Vdist, Normal_cv, 6.137241, 0.2);

Simulation { # 1:
  Print(Dose, 1);
}
END.
