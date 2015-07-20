#!/usr/bin/env gnuplot

# Use the parula color palette
set palette defined (0 '#352a87', 1 '#0363e1', 2 '#1485d4', 3 '#06a7c6', 4 '#38b99e', 5 '#92bf73', 6 '#d9ba56', 7 '#fcce2e', 8 '#f9fb0e')

# Define plot options
set pm3d map
set pm3d interpolate 0,0

set xlabel  'Energy ϵ/Δ'   offset 0,0
set ylabel  'Position z/ξ' offset 0,0

set xrange  [-1.5:1.5]
set cbrange [ 0.0:2.0]

set xtics  -1.5,0.5,1.5
set cbtics -10,1,10
set ytics  -10,1,10

set format x '%3.1f'

# Plot a file provided from stdin to interactive terminal
set terminal wxt enhanced font 'Verdana,10'
splot '/dev/stdin' using 2:1:3 notitle

# Wait for a mouse click or keyboard press, and loop
pause mouse

# Replot the same file to the file ~/gnuplot.png
set terminal pngcairo size 12000,10000 crop enhanced font 'Helvetica,250'
set output '~/gnuplot.png'
splot '/dev/stdin' using 2:1:3 notitle

# Loop
reread

