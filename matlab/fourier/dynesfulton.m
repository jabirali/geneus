function dynesfulton(inputfile, threshold)
  % This function uses the Dynes-Fulton method [PRB 3(9) 3015 (1971)] to estimate
  % the current density distribution based on critical current measurements.

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % IMPORT AND PREPROCESS THE CURRENT-FIELD MEASUREMENTS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Import data from an input file
  data    = load(inputfile);
  field   = data(:,1);
  current = data(:,2);

  % Remove any signs and direct currents
  current = abs(current);
  current = current - min(current);

  % Estimate the local curvature
  c = curvature(current, field);

  % Identify spikes in the curvature
  u = trim(spikes(c, threshold));

  % Create a sign array that flips after each spike
  s = (-1) .^ cumsum(u);

  % Use this to restore the sign of the measured current
  current = s .* current;

  % Make sure the central lobe is a peak
  if (abs(max(current)) < abs(min(current)))
    current = -current;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % INTERPOLATE, FIND THE CURRENT DENSITY, AND SMOOTH
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Interpolate the input data to a higher precision
  field_i   = linspace(min(field), max(field), 10240);
  current_i = pchip(field, current, field_i);

  % Perform a Fourier transform of the results
  [density, position] = fourier(current, field);

  % Smooth the resulting current density to get the trend
  density_s = smooth(position, density, 'lowess');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % VISUALIZE THE FINAL RESULTS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % First subfigure
  figure;
  subplot(3,1,1);

  % Plot the original measurements
  plot(data(:,1), data(:,2), 'b.', data(:,1), data(:,2), 'k-');
  xlabel('H');
  ylabel('J_c(H)');

  % Second subfigure
  subplot(3,1,2);

  % Plot the reconstructed current
  current_p = current; current_p(current < 0) = NaN;
  current_m = current; current_m(current > 0) = NaN;
  plot(field_i, current_i, 'k-', field, current_p, 'b.', field, current_m, 'r.');
  xlabel('H');
  ylabel('J(H)');

  % Third subfigure
  subplot(3,1,3);

  % Plot the Fourier transformation
  hold on
  area(position, density,   'FaceColor', 'k');
  plot(position, density_s, 'r-');
  hold off
  xlabel('x');
  ylabel('j(x)');
end

function v = trim(u)
  % Takes a boolean array, and replaces consecutive true-valued regions in the array
  % with a single true element in the center of the region. For instance:
  %  trim([0,1,1,1,0,0,0,1,1,1,1,1,0,0]) = [0,0,1,0,0,0,0,0,0,1,0,0,0,0]

  % Differentiate the boolean array — the result should be +1 at the start of a
  % true-valued region, -1 at the end of such a region, and 0 otherwise
  d = diff(u);

  % Identify the start and stop of true-valued regions
  p = find(d > 0);
  q = find(d < 0);

  % Finally, construct the trimmed boolean array
  v = zeros(size(u));
  for n=1:min(size(p),size(q))
    v(ceil(1+(p(n)+q(n))/2)) = 1;
  end
end

function u = spikes(y, t)
  % Identifies maxima in a dataset that can be classified as outliers.
  % [Assumption: the original data set a median that is roughly zero.]

  u = y > t*mad(y,1);
end

function c = curvature(y, x)
  % Calculates the local curvature of some function y(x), using finite-difference
  % approximations for the numerical derivatives dy/dx and d²y/dx² of the arrays.

  % Calculate the first and second derivatives
  d1 = differentiate1(y, x);
  d2 = differentiate2(y, x);

  % Calculate the signed local curvature 
  c = d2 ./ (1+d1.^2).^(3/2);
end

function d = differentiate1(y, x)
  % Calculates the first-derivative dy/dx, using a central-difference approximation in 
  % the interior domain, and a forward/backward-difference approximation at the edges.

  % Initialize variables
  N = min(length(x),length(y));
  d = zeros(N,1);

  % Calculate the finite-difference first-derivative
  d(  1  ) = (y( 2 )-y(  1  )) ./ (x( 2 )-x(  1  ));
  d(2:N-1) = (y(3:N)-y(1:N-2)) ./ (x(3:N)-x(1:N-2));
  d(  N  ) = (y( N )-y( N-1 )) ./ (x( N )-x( N-1 ));
end

function d = differentiate2(y, x)
  % Calculates the second-derivative d²y/dx², using a central-difference approximation
  % in the interior domain of the function y(x), and skipping calculation at the edges.

  % Initialize variables
  N = min(length(x),length(y));
  d = zeros(N,1);

  % Calculate the finite-difference second-derivative
  d(2:N-1) = 4 * (y(3:N) - 2*y(2:N-1) + y(1:N-2)) ./ (x(3:N) - x(1:N-2)).^2;
end

function [Y,X] = fourier(y, x)
  % Performs a Fourier transformation of an array y(x).

  % Perform the transformation
  Y = fft(y);

  % Reorganize the results
  p = floor(1+length(y)/2);
  Y = abs(real([Y(p+1:end); Y(1:p)]));

  % Construct a fitting domain
  q = 2/(mean(diff(x)));
  X = linspace(-q, q, length(Y))';

  % Calculate the norm using trapezoid integration
  s = sqrt(sum( (Y(2:end).^2+Y(1:end-1).^2) .* (X(2:end)-X(1:end-1)))/2);

  % Normalize the output result
  Y = Y/s;
end
