function dynesfulton(inputfile, outputfile)
  % This function uses the Dynes-Fulton method [PRB 3(9) 3015 (1971)] to estimate
  % the current density distribution based on critical current measurements.

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % IMPORT AND PREPROCESS THE CURRENT-FIELD MEASUREMENTS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Import data from an input file
  data    = load(inputfile);

  % Extract data using SI units
  field   = data(:,1);
  current = data(:,2);



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % INTERPOLATE, FIND THE CURRENT DENSITY, AND SMOOTH
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Remove duplicates
  d = diff(field);
  field_u   = [field(1)  ];
  current_u = [current(1)];
  for n=1:size(d)
	  if (d(n) > 0)
		  field_u   = [field_u;   field(n+1)  ];
		  current_u = [current_u; current(n+1)];
	  end
  end

  % Interpolate the input data to a higher precision
  field_i   = linspace(-300, +300, 10240)';
  current_i = pchip(field_u, current_u, field_i);
  current_j = (current_i + fliplr(current_i')')/2;
  current_k = [current_i(1:5120); fliplr(current_i(1:5120)')'];

  % Perform a Fourier transform of the results
  beta_i                  = field_i * 0.48359785;
  [density_i, position_i] = fourier(current_i, beta_i);
  [density_j, position_j] = fourier(current_j, beta_i);
  [density_k, position_k] = fourier(current_k, beta_i);

  % Smooth the resulting current density to get the trend
  %density = smooth(position, density, 'lowess');



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % VISUALIZE THE FINAL RESULTS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Plot the reconstructed current
  subplot(3,2,1);
  plot(field_i, current_i, 'k-');
  title('Asymmetric data');
  xlabel('Applied field μ₀H [mT]');
  ylabel('Total current J(H)');
  xlim([-300, +300]);
  ylim([-1.1, +1.1]);

  subplot(3,2,3);
  plot(field_i, current_j, 'k-');
  title('Symmetrized data');
  xlabel('Applied field μ₀H [mT]');
  ylabel('Total current J(H)');
  xlim([-300, +300]);
  ylim([-1.1, +1.1]);

  subplot(3,2,5);
  plot(field_i, current_k, 'k-');
  title('Mirrored data');
  xlabel('Applied field μ₀H [mT]');
  ylabel('Total current J(H)');
  xlim([-300, +300]);
  ylim([-1.1, +1.1]);

  % Plot the Fourier transformation
  subplot(3,2,2);
  area(position_i, density_i,   'FaceColor', 'k');
  title('Asymmetric data');
  xlabel('Position x [μm]');
  ylabel('Current density j(x)');
  xlim([-0.5, +0.5]);

  subplot(3,2,4);
  area(position_j, density_j,   'FaceColor', 'k');
  title('Symmetrized data');
  xlabel('Position x [μm]');
  ylabel('Current density j(x)');
  xlim([-0.5, +0.5]);

  subplot(3,2,6);
  area(position_k, density_k,   'FaceColor', 'k');
  title('Mirrored data');
  xlabel('Position x [μm]');
  ylabel('Current density j(x)');
  xlim([-0.5, +0.5]);

  % Save the results
  output = [position_i, density_i];
  save(outputfile, 'output', '-ascii');
end

function [Y,X] = fourier(y, x)
  % Performs a Fourier transformation of an array y(x).

  % Perform the transformation
  Y = abs(fftshift(ifft(y)));

  % Construct a fitting domain
  q = 1/(mean(diff(x)));
  X = linspace(-q/2, +q/2, length(Y))';

  % Calculate the norm using trapezoid integration
  s = sqrt(sum( (Y(2:end).^2+Y(1:end-1).^2) .* (X(2:end)-X(1:end-1)))/2);

  % Normalize the output result
  Y = Y/s;
end
