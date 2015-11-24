#!/usr/bin/env gnuplot

# Parula color palette
set palette defined (0 '#352a87', 1 '#0363e1', 2 '#1485d4', 3 '#06a7c6', 4 '#38b99e', 5 '#92bf73', 6 '#d9ba56', 7 '#fcce2e', 8 '#f9fb0e')

# Configure the borders
set style line 11 lc rgb '#808080' lt 1
set tics nomirror out scale 0.75

# Configure the axes
set xlabel  'Interface polarization P'        offset 0,0
set ylabel  'Spin-mixing conductance G_φ/G_0' offset 0,0

set xrange  [-1.001:1.001]
set yrange  [-2.001:2.001]
set cbrange [-0.001:2.001]

set xtics  -1.0,0.5,1.0
set ytics  -2.0,1.0,2.0
set cbtics -10,1,10

set format x '%3.1f'
set format y '%3.1f'

# Plot a file provided from stdin to interactive terminal
#set terminal wxt enhanced font 'Gillius ADF,12'
set terminal pdf enhanced font 'Gillius ADF,16' linewidth 1.5
set out 'params_spinactive.pdf'
plot '/dev/stdin' using ($1/200-0.99575):($2/100-2):3 matrix notitle with image

# Wait for a mouse click
#pause mouse
