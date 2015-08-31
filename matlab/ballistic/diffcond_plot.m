function diffcond_plot(transmission, polarization, scattering, temperature)
  % Calculates and plots the differential conductance as function of voltage for
  % a ballistic superconductor/normal-metal bilayer with a spin-active tunneling
  % interface. The temperature, transmission probability, spin-polarization, and
  % inelastic scattering rate must be specified as arguments. The temperature is
  % normalized to the superconducting critical temperature Tc, and the inelastic
  % scattering is an imaginary energy term normalized to the superconducting gap.



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                      PREPARATIONS FOR THE CALCULATIONS                     %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Parameters to use in calculations
  voltage = linspace(-3,+3,200);
  energy  = linspace(-6,+6,400);
  angle   = linspace(0,pi/2,1000);

  % Spin-dependent transmission probabilities
  D_up = transmission * (1+polarization);
  D_dn = transmission * (1-polarization);

  % Spin-dependent reflection probabilities
  R_up = 1-D_up;
  R_dn = 1-D_dn;

  % Coefficients required in the equations below
  U = R_up*R_dn;
  V = -2*sqrt(U);
  A = 2*D_up*D_dn/(D_up+D_dn);
  B = (1+U)/A;
  C = V/A;



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                          SPECTRAL CURRENT DENSITY                          %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  spectral = zeros(length(energy),length(temperature));
  for k=1:length(temperature)
    % Calculate the superconducting gap at this temperature
    gap = tanh(1.74*sqrt(1/(temperature(k)+1e-16)-1));

    % Calculate the spectral current density as a function of energy
    for a = 1:length(angle)
      phaseshift  = pi * exp(-4*sin(angle(a))^2);
      probability = cos(angle(a)) * (range(angle)/length(angle));
      for m=1:length(energy)
        if (abs(energy(m)) < gap)
          % Calculation below the superconducting gap
          % [eqs. (45) and (50) in PRB 70 134510]
          delta         = acos((energy(m)+i*scattering)/gap);
          spectral(m,k) = spectral(m,k)                       ...
                        + probability                         ...
                        * ( 1/(B+C*cos(2*delta+pi*phaseshift))...
                          + 1/(B+C*cos(2*delta-pi*phaseshift)) );
        else
          % Calculation above the superconducting gap
          % [eqs. (46) and (49) in PRB 70 134510]
          delta         = acosh((abs(energy(m))+i*scattering)/gap);
          spectral(m,k) = spectral(m,k)                                            ...
                        + probability                                              ...
                        * ( 2*cosh(delta)*(A*exp(-delta) + 2*sinh(delta))          ...
                          / (exp(2*delta) + U*exp(-2*delta) + V*cos(pi*phaseshift)) );
        end
      end
    end
  end



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                            TOTAL CURRENT DENSITY                           %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  current  = zeros(length(voltage),length(temperature));
  step     = range(energy)/(length(energy)-1);
  for k=1:length(temperature)
    for n=1:length(voltage)
      for m=1:length(energy)
        % Calculate the total current density as a function of voltage bias
        % [eq.(51) in PRB 70 134510, where 0.8819... is half the BCS ratio]
        current(n,k) = current(n,k)                                                            ...
                     + real(spectral(m,k)) * step                                              ...
                     * (tanh(0.8819384944310228*(energy(m)+voltage(n))/(temperature(k)+1e-16)) ...
                       -tanh(0.8819384944310228*(energy(m)-voltage(n))/(temperature(k)+1e-16)))/4;
      end 
    end
  end



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                            DIFFERENTIAL CONDUCTANCE                        %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Calculate the differential conductance as a function of voltage
  diffcond = zeros(length(voltage)-1,length(temperature));
  for k=1:length(temperature)
    diffcond(:,k) = diff(current(:,k)) ./ diff(voltage(:));
  end

  % Plot the differential conductance (overlay)
  figure;
  hold on;
  plot((voltage(1:end-1)+voltage(2:end))/2, diffcond);
  xlabel('Voltage eV/\Delta_0');
  ylabel('Differential conductance R·dI/dV');
  legend(strtrim(cellstr(strcat('T=',num2str(temperature')))))
  axis([min(voltage) max(voltage) -2 4]);

  % Plot the differential conductance (gradient)
  figure;
  hold on;
  surf((voltage(1:end-1)+voltage(2:end))/2, temperature, diffcond', 'EdgeColor', 'None');
  hcb = colorbar;
  xlabel('Voltage eV/\Delta_0');
  ylabel('Temperature T/T_c');
  caxis([0 2]);
  set(hcb,'YTick',[0,1,2]);
  colormap(parula(256));
end
