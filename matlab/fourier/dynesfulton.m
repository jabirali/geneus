function dynesfulton(inputfile, outputfile)
  % Load input file
  data    = load(inputfile);
  field   = data(:,1);
  current = data(:,2);
  N       = length(field);

  % Estimate the first derivative
  d1 = zeros(N,1);
  for n=(1+1):(N-1)
    d1(n) = (current(n+1)-current(n-1))/(field(n+1)-field(n-1));
  end
  d1(1) = (current(1+1)-current(1))/(field(1+1)-field(1));
  d1(N) = (current(N)-current(N-1))/(field(N)-field(N-1));

  % Create the indicator function
  % TODO: Make this more robust using med+mad to find d1 spikes
  u = ones(N,1);
  for n=2:N
    if (d1(n) > 0 && d1(n-1) < 0)
      u(n:N) = -u(n:N);
    end
  end

  % Multiply by the indicator
  current = current .* u;

  % Make sure the central peak is a maximum
  if (abs(max(current)) < abs(min(current)))
    current = -current;
  end

  % Perform a Fourier transform
  % TODO: Make this part work
  fourier = ifft(field)

  % Dump processed output
  data = [field, current];
  save(outputfile, 'data', '-ascii');
end
