% Grid resolution
points = 256;

% Normalized magnetic field
limit = pi/2 + 20*pi;
field = linspace(-limit, +limit, points)';

% Fraunhofer pattern
current = zeros(points,1);
for n=1:points
  current(n) = (abs(sin(field(n)/2))+1e-16)/(abs(field(n)/2)+1e-16);
end 

% Dump output
data = [field, current];
save fraunhofer.dat data -ascii
