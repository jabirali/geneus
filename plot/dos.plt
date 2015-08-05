#!/usr/bin/env gnuplot

# Parula color palette
set palette defined (0 '#352a87', 1 '#0363e1', 2 '#1485d4', 3 '#06a7c6', 4 '#38b99e', 5 '#92bf73', 6 '#d9ba56', 7 '#fcce2e', 8 '#f9fb0e')

# Configure the borders
set style line 11 lc rgb '#808080' lt 1
set tics nomirror out scale 0.75

# Configure the axes
set xlabel  'Energy ϵ/Δ'   offset 0,0
set ylabel  'Position z/ξ' offset 0,0

set xrange  [-1.501:1.501]
set cbrange [ 0.0:2.0]

set xtics  -1.501,0.5,1.501
set cbtics -10,1,10
set ytics  -10,0.49999,10

set format x '%3.1f'
set format y '%3.1f'

# Configure the heatmap
set pm3d map
set pm3d interpolate 0,0

# Plot a file provided from stdin to interactive terminal
set terminal wxt enhanced font 'Gillius ADF,12'
splot '/dev/stdin' using 2:1:3 notitle

# Wait for a mouse click
pause mouse

# Rerun the file
reread
