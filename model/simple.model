States = {y0, y1, y2};

Outputs = {yout};

# Parameters
k1 = 1;
k2 = 1;
k3 = 1;

Dynamics {
  
  dt(y0) = -k1 * y0 + k2 * y1 * y2; 
  dt(y2) =  k3 * y1 * y1; 
  dt(y1) = -dt(y0) - dt(y2);
  
  yout = y0 + y1 + y2;
  
} # End of Dynamics

Events {
  y0 = 1;
}

Roots { # here we need inlining, otherwise 'gout' will not be understood 
  Inline ( gout[0] = y[0] - 0.5; );
}

End.