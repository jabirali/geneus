#!/usr/bin/env gnuplot

# Configure the borders
set style line 11 lc rgb '#808080' lt 1
set tics nomirror out scale 0.75

# Configure the axes
set xtics  0.00,0.25,1.00
set xrange [0.00:1.00]

set format x '%3.2f'

# Plot the electric current
set terminal wxt
set xlabel 'Phase difference φ/π'
set ylabel 'Critical current I_c/I_c_0'
plot filename u 1:2 notitle ls 99

# Wait for a mouse click
pause mouse
