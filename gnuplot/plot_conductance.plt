#!/usr/bin/env gnuplot

# Configure the borders
set style line 11 lc rgb '#808080' lt 1
set tics nomirror out scale 0.75

# Configure the axes
set xlabel 'Voltage eV/Δ'                   offset 0,0
set ylabel 'Differential conductance dI/dV' offset 0,0

set xtics -1.501,0.5,1.501
set ytics -10,1,10

set xrange [-1.501:1.501]
#set yrange [ 0.0:4.0]

set format x '%3.1f'
set format y '%3.1f'

# Parula colormap
set style line 11 lt 1 lc rgb '#352a87' # blue
set style line 12 lt 1 lc rgb '#d95319' # orange
set style line 13 lt 1 lc rgb '#edb120' # yellow
set style line 14 lt 1 lc rgb '#7e2f8e' # purple
set style line 15 lt 1 lc rgb '#77ac30' # green
set style line 16 lt 1 lc rgb '#4dbeee' # light-blue
set style line 17 lt 1 lc rgb '#a2142f' # red

# Plot a file provided from stdin to interactive terminal
set style fill transparent solid 0.5
set style data lines
set terminal wxt enhanced font 'Gillius ADF,12'
plot '/dev/stdin' notitle ls 11 lw 2 with filledcurves x1

# Wait for a mouse click
pause mouse
