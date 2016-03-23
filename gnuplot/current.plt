#!/usr/bin/env gnuplot

# Configure the borders
set style line 11 lc rgb '#808080' lt 1
set tics nomirror out scale 0.75

# Configure the axes
set xlabel 'Phase difference φ/π' offset 0,0
set ylabel 'Current'              offset 0,0

set xtics  0.00,0.25,1.00
set xrange [0.00:1.00]

set format x '%3.2f'

# Plot a file provided from stdin to interactive terminal
set terminal wxt enhanced font 'Gillius ADF,12'
plot filename u 1:2 title 'I_e' ls 99, \
     filename u 1:3 title 'I_{/Symbol s}_,_x' ls  1, \
     filename u 1:4 title 'I_{/Symbol s}_,_y' ls  3, \
     filename u 1:5 title 'I_{/Symbol s}_,_z' ls  5

# Wait for a mouse click
pause mouse
