function c = dynesfulton(inputfile, outputfile)
  % Load input file
  data    = load(inputfile);
  field   = data(:,1);
  current = abs(data(:,2));

  % Estimate the local curvature
  c = curvature(current, field);

  % Identify spikes in the curvature
  u = spikes(c);

  % Create a sign array that flips after each spike
  s = (-1) .^ cumsum(u);

  % Use this to restore the sign of the measured current
  current = s .* current;

  % Make sure the central peak is a maximum
  if (abs(max(current)) < abs(min(current)))
    current = -current;
  end

  % Perform a Fourier transform
  % TODO: Make this part work
  fourier = ifft(field)

  % Dump processed output
  %data = [field, u];
  data = u;
  save(outputfile, 'data', '-ascii');
end

function u = spikes(y)
  % Uses an outlier-detection algorithm to identify abberant maxima in a dataset.

  % Find the median and median-absolute-deviation of the input array
  N = length(y);
  m = median(y);
  s = mad(y,1);

  % Identify regions that exceed a statistical threshold value
  u = y > m + 3*s;

  % Check where the outlier regions start and stop
  %d  = differentiate1(u, (1:N)');
  %u1 = find(d > 0);
  %u2 = find(d < 0);

  %um = mean([u1, u2]')';

  %um = u;
end

function c = curvature(y, x)
  % Calculates the local curvature of some function y(x), using a finite-difference
  % approximation for the numerical derivatives dy/dx and d²y/dx² of the arrays.

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
  d(2:N-1) = (y(3:N) - 2*y(2:N-1) + y(1:N-2)) ./ (x(3:N) - x(1:N-2)).^2;
end
