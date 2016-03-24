#!/usr/bin/env gnuplot

# Plot the gap magnitude
set terminal wxt 
set xlabel 'Position z/ξ'
set ylabel 'Superconducting gap Δ/Δ_0'
plot filename u 1:2 notitle ls 99

# Wait for a mouse click
pause mouse

# Plot the gap phase
set terminal wxt 
set xlabel 'Position z/ξ'
set ylabel 'Superconducting phase φ/π'
plot filename u 1:3 notitle ls 99

# Wait for a mouse click
pause mouse
