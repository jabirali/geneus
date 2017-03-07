% Grid resolution
points = 256;

% Normalized magnetic field
field = (pi/2) * linspace(-40, +40, points)';

% Fraunhofer pattern
current = zeros(points,1);
for n=1:points
  current(n) = 1 + abs(sin(field(n)/2))/(abs(field(n)/2)+1e-16);
end 

% Dump output
data = [field, current];
save fraunhofer.dat data -ascii
