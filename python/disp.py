#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys               import argv, exit
from numpy             import *
from matplotlib.pyplot import *

# Define color palette and colormap
palette = {'heatmap' : 'inferno',
           'black'   : '#000000',
           'blue'    : '#377eb8',
           'green'   : '#4daf4a',
           'red'     : '#e41a1c'}

def plot_density(data):
    """Function for plotting the density of states."""

    # Extract and reshape data
    position = unique(data[:,0])
    energy   = unique(data[:,1])
    density  = data[:,2].reshape(len(position), len(energy))

    # Plot the densirt of states
    figure()
    gcf().canvas.set_window_title('Density of states')
    pcolormesh(energy, position, density, vmin = 0.0, vmax = 2.0, cmap = palette['heatmap'])
    colorbar(ticks = [0.0,0.5,1.0,1.5,2.0])
    axis([-2.0, +2.0, position.min(), position.max()])
    xlabel(r'Energy $\epsilon/\Delta$')
    ylabel(r'Position $z/\xi$')
    show()

def plot_current(data):
    """Function for plotting the charge and spin currents."""

    # Extract and reshape data
    position = unique(data[:,0])
    charge   = data[:,1]
    spin_x   = data[:,2]
    spin_y   = data[:,3]
    spin_z   = data[:,4]

    # Plot the charge current
    figure()
    gcf().canvas.set_window_title('Charge current')
    plot(position, charge, palette['black'], linewidth = 2)
    axis([position.min(), position.max(), min(0.0,charge.min())*1.05, max(0.0,charge.max())*1.05])
    ylabel(r'Charge current $I_{e}/I_{e0}$')
    xlabel(r'Position $z/\xi$')
    show()

    # Plot the spin current
    figure()
    gcf().canvas.set_window_title('Spin current')
    plot(position, spin_x, palette['blue'], 
         position, spin_y, palette['green'], 
         position, spin_z, palette['red'], 
         linewidth = 2)
    axis([position.min(), 
          position.max(), 
          min(0.0, spin_x.min(), spin_y.min(), spin_z.min())*1.05, 
          max(0.0, spin_x.max(), spin_y.max(), spin_z.max())*1.05])
    ylabel(r'Spin current $I_{\sigma}/I_{\sigma 0}$')
    xlabel(r'Position $z/\xi$')
    legend(['Spin-$x$', 'Spin-$y$', 'Spin-$z$'])
    show()

def plot_gap(data):
    """Function for plotting the superconducting gap and phase."""

    # Extract and reshape data
    position = unique(data[:,0])
    gap      = data[:,1]
    phase    = data[:,2]

    # Plot the superconducting gap
    figure()
    gcf().canvas.set_window_title('Superconducting gap')
    plot(position, gap, palette['black'], linewidth = 2)
    axis([position.min(), position.max(), min(-1e-6,gap.min())*1.05, max(+1e-6,gap.max())*1.05])
    ylabel(r'Superconducting gap $\Delta/\Delta_0$')
    xlabel(r'Position $z/\xi$')
    show()

    # Plot the superconducting phase
    figure()
    gcf().canvas.set_window_title('Superconducting phase')
    plot(position, phase, palette['black'], linewidth = 2)
    axis([position.min(), position.max(), min(-1e-6,phase.min())*1.05, max(+1e-6,phase.max())*1.05])
    ylabel(r'Superconducting phase $\varphi/\pi$')
    xlabel(r'Position $z/\xi$')
    show()

if __name__ == "__main__":
    # Check command line arguments
    try:
        filename = argv[1]
    except:
        print('Usage: {} <filename.dat>'.format(argv[0]))
        exit(0)

    # Read data from file
    with open(filename, 'r') as f:
        desc = f.readline()   # File description
        data = loadtxt(f)     # File contents

    # Split the file description into fields
    field = [q.strip() for q in desc[2:].split('\t')]

    # Call the appropriate plotting method
    if field == ['Position', 'Energy', 'Density of states']:
        plot_density(data)
    elif field == ['Position', 'Charge current', 'Spin-x current', 'Spin-y current', 'Spin-z current']:
        plot_current(data)
    elif field == ['Position', 'Gap magnitude', 'Gap phase']:
        plot_gap(data)
    else:
        print('Sorry, unrecognized data format...')
