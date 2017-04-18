function dynesfulton(inputfile, threshold)
  % This function uses the Dynes-Fulton method [PRB 3(9) 3015 (1971)] to estimate
  % the current density distribution based on critical current measurements.

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % IMPORT AND PREPROCESS THE CURRENT-FIELD MEASUREMENTS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Import data from an input file
  data    = load(inputfile);
  field   = data(:,1);
  current = abs(data(:,2));

  % Identify nodepoints in the current
  u = trim(current < threshold*max(current))
  % TODO: If this doesn't work well, try calculating the local
  %       median m and median-absolute-deviation s, where local
  %       here means calculated from the e.g. 20 neighbouring
  %       elements. If the point is below m-3s, then we have
  %       a local minimum, and can use trim to find its center.

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
  %field   = linspace(min(field), max(field), 10240)';
  %current = pchip(field, current, field_i);

  % Perform a Fourier transform of the results
  [density, position] = fourier(current, field);

  % Smooth the resulting current density to get the trend
  %density = smooth(position, density, 'lowess');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % VISUALIZE THE FINAL RESULTS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % First subfigure
  figure;
  subplot(3,1,1);

  % Plot the original measurements
  plot(data(:,1), data(:,2), 'k.-');
  xlabel('H');
  ylabel('J_c(H)');

  % Second subfigure
  subplot(3,1,2);

  % Plot the reconstructed current
  current_p = current; current_p(current < 0) = NaN;
  current_m = current; current_m(current > 0) = NaN;
  plot(field, current, 'k-', field, current_p, 'b.', field, current_m, 'r.');
  xlabel('H');
  ylabel('J(H)');

  % Third subfigure
  subplot(3,1,3);

  % Plot the Fourier transformation
  area(position, density,   'FaceColor', 'k');
  xlabel('x');
  ylabel('j(x)');
  xlim([-3, 3]);
end

function v = trim(u)
  % Takes a boolean array, and replaces consecutive true-valued regions in the 
  % array with a single true element in the center of the region. For instance:
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
