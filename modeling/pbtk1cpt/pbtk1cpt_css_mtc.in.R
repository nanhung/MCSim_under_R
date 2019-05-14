MonteCarlo ("", 1000, 10101010);

Distrib (Fgutabs, Uniform, 0.8, 1.0);
Distrib (kelim, Normal_cv, 0.02283233, 0.2);
Distrib (Vdist, Normal_cv, 6.137241, 0.2);

Simulation { # 1:
  Print(css, 0.228291);
}
END.
