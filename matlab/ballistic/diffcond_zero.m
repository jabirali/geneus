function result = diffcond_zero(voltage, transmission, polarization, phaseshift)
  % Calculates the differential conductance as a function of the normalized
  % voltage eV/\Delta for a ballistic superconductor/normal-metal bilayer with a
  % spin-active interface. The transmission probability, spin-polarization, and 
  % spin-dependent phase shifts at the spin-active interface must be specified.
  
  % Calculate the spin-dependent transmission probabilities
  D_up = transmission * (1+polarization);
  D_dn = transmission * (1-polarization);

  % Calculate the spin-dependent reflection probabilities
  R_up = 1-D_up;
  R_dn = 1-D_dn;

  % Calculate some auxiliary quantities used in the equations
  A =  2*D_up*D_dn/(D_up+D_dn);
  B = (1 + R_up*R_dn)/A;
  C = (-2*sqrt(R_up*R_dn))/A;
  U = R_up*R_dn;
  V = (-2*sqrt(R_up*R_dn));
  
  % Calculate the differential conductance
  for n=1:length(voltage)
    if (abs(voltage(n)) < 1)
      % Below the superconducting gap
      delta     = acos(voltage(n));
      result(n) = 1/(B+C*cos(2*delta+phaseshift))                            ...
                + 1/(B+C*cos(2*delta-phaseshift));
    else
      % Above the superconducting gap
      delta     = acosh(abs(voltage(n)));
      result(n) = 2*cosh(delta)*(A*exp(-delta) + 2*sinh(delta))              ...
                / (exp(2*delta) + U*exp(-2*delta) + V*cos(phaseshift));
    end
  end
end
