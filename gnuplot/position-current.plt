#!/usr/bin/env gnuplot

# Configure the borders
set style line 11 lc rgb '#808080' lt 1
set tics nomirror out scale 0.75

# Plot the electric current
set terminal wxt
set xlabel 'Position z/ξ'
set ylabel 'Charge current I_e/I_e_0'
plot filename u 1:2 notitle ls 99

# Wait for a mouse click
pause mouse

# Plot the spin current
set terminal wxt
set xlabel 'Position z/ξ'
set ylabel 'Spin current I_s/I_s_0'
plot filename u 1:3 title 'I_x' ls 3, \
     filename u 1:4 title 'I_y' ls 5, \
     filename u 1:5 title 'I_z' ls 1

# Wait for a mouse click
pause mouse
