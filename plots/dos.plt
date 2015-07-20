#!/usr/bin/env gnuplot

set terminal wxt enhanced font 'Verdana,10'

#@png
#set output "output.png"
#plot x ls 1, -x ls 2, x**3 ls 3

#print GPVAL_TERMINALS
#set terminal png size 4096,4096 crop enhanced font "/usr/share/fonts/truetype/times.ttf,100" dashlength 2

#set terminal tikz standalone
#set output 'output.tex'

#set terminal svg size 350,262 fname 'Verdana' fsize 10
#set output 'output.svg'

#set terminal pngcairo size 10000,10000 crop enhanced font 'Helvetica' fontscale 1.0
#set output 'output.png'

# Use the parula color palette
set palette defined (0 '#352a87', 1 '#0363e1', 2 '#1485d4', 3 '#06a7c6', 4 '#38b99e', 5 '#92bf73', 6 '#d9ba56', 7 '#fcce2e', 8 '#f9fb0e')


set pm3d map
set pm3d interpolate 0,0

set xlabel 'Energy ϵ/Δ'   offset 0,2
set ylabel 'Position z/ξ' offset 0

set cbrange [ 0.0:2.0]
set xrange  [-1.5:1.5]
#set yrange  auto

set cbtics 0,1,2 #offset -0.70,0.15
unset xtics
unset ytics

splot '/dev/stdin' using 2:1:3 notitle

pause mouse
reread
#replot
