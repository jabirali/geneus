% Grid resolution
points = 128;

% Normalized magnetic field
field = (pi/2) * linspace(-10, +10, points)';

% Fraunhofer pattern
current = zeros(points,1);
for n=1:points
  current(n) = abs(sin(field(n)))/(abs(field(n))+1e-16);
end 

% Dump output
data = [field, current];
save fraunhofer.dat data -ascii
