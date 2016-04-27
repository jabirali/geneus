#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys               as s
import numpy             as n
import matplotlib.pyplot as p

# Define plot style
p.style.use('ggplot')

# Define color palette and colormap
palette = {'heatmap' : 'viridis',
           'black'   : '#000000',
           'blue'    : '#377eb8',
           'green'   : '#4daf4a',
           'red'     : '#e41a1c'}

# Define font type and size
fontset = {'fontname' : 'Sans',
           'size'     : '16'  }

def plot_density(data):
    """Function for plotting the density of states."""

    # Extract and reshape data
    position = n.unique(data[:,0])
    energy   = n.unique(data[:,1])
    density  = data[:,2].reshape(len(position), len(energy))

    # Plot the density of states
    p.figure()
    p.gcf().canvas.set_window_title('Density of states')
    p.pcolormesh(energy, position, density, vmin = 0.0, vmax = 2.0, cmap = palette['heatmap'])
    p.colorbar(ticks = [0.0,0.5,1.0,1.5,2.0])
    p.axis([-2.0, +2.0, position.min(), position.max()])
    p.xlabel(r'Energy $\epsilon/\Delta$', **fontset)
    p.ylabel(r'Position $z/\xi$', **fontset)
    p.show()

def plot_current(data):
    """Function for plotting the charge and spin currents."""

    # Extract and reshape data
    position = n.unique(data[:,0])
    charge   = data[:,1]
    spin_x   = data[:,2]
    spin_y   = data[:,3]
    spin_z   = data[:,4]

    # Plot the charge current
    p.figure()
    p.gcf().canvas.set_window_title('Charge current')
    p.plot(position, charge, palette['black'], linewidth = 2)
    p.axis([position.min(), position.max(), min(0.0,charge.min())*1.05, max(0.0,charge.max())*1.05])
    p.ylabel(r'Charge current $I_{e}/I_{e0}$', **fontset)
    p.xlabel(r'Position $z/\xi$', **fontset)
    p.show()

    # Plot the spin current
    p.figure()
    p.gcf().canvas.set_window_title('Spin current')
    p.plot(position, spin_x, palette['blue'], 
           position, spin_y, palette['green'], 
           position, spin_z, palette['red'], 
           linewidth = 2)
    p.axis([position.min(), 
            position.max(), 
            min(0.0, spin_x.min(), spin_y.min(), spin_z.min())*1.05, 
            max(0.0, spin_x.max(), spin_y.max(), spin_z.max())*1.05])
    p.ylabel(r'Spin current $I_{\sigma}/I_{\sigma 0}$', **fontset)
    p.xlabel(r'Position $z/\xi$', **fontset)
    p.legend(['Spin-$x$', 'Spin-$y$', 'Spin-$z$'])
    p.show()

def plot_gap(data):
    """Function for plotting the superconducting gap and phase."""

    # Extract and reshape data
    position = n.unique(data[:,0])
    gap      = data[:,1]
    phase    = data[:,2]

    # Plot the superconducting gap
    p.figure()
    p.gcf().canvas.set_window_title('Superconducting gap')
    p.plot(position, gap, palette['black'], linewidth = 2)
    p.axis([position.min(), position.max(), min(-1e-6,gap.min())*1.05, max(+1e-6,gap.max())*1.05])
    p.ylabel(r'Superconducting gap $\Delta/\Delta_0$', **fontset)
    p.xlabel(r'Position $z/\xi$', **fontset)
    p.show()

    # Plot the superconducting phase
    p.figure()
    p.gcf().canvas.set_window_title('Superconducting phase')
    p.plot(position, phase, palette['black'], linewidth = 2)
    p.axis([position.min(), position.max(), min(-1e-6,phase.min())*1.05, max(+1e-6,phase.max())*1.05])
    p.ylabel(r'Superconducting phase $\varphi/\pi$', **fontset)
    p.xlabel(r'Position $z/\xi$', **fontset)
    p.show()

if __name__ == "__main__":
    # Check command line arguments
    try:
        filename = s.argv[1]
    except:
        print('Usage: {} <filename.dat>'.format(s.argv[0]))
        s.exit(0)

    # Read data from file
    with open(filename, 'r') as f:
        desc = f.readline()   # File description
        data = n.loadtxt(f)   # File contents

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
