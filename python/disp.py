#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Define plot style
plt.style.use('ggplot')
mpl.rcParams['figure.figsize']   = (10, 6)
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family']      = 'STIXGeneral'
mpl.rcParams['font.size']        = '16'
plt.rcParams['axes.labelsize']    = plt.rcParams['font.size']
plt.rcParams['axes.titlesize']    = plt.rcParams['font.size'] * 1.5
plt.rcParams['legend.fontsize']   = plt.rcParams['font.size']
plt.rcParams['xtick.labelsize']   = plt.rcParams['font.size']
plt.rcParams['ytick.labelsize']   = plt.rcParams['font.size']
plt.rcParams['savefig.dpi']       = 600
plt.rcParams['xtick.major.size']  = 0
plt.rcParams['xtick.minor.size']  = 0
plt.rcParams['ytick.major.size']  = 0
plt.rcParams['ytick.minor.size']  = 0
plt.rcParams['legend.frameon']    = False
plt.rcParams['legend.loc']        = 'upper center'
plt.rcParams['legend.borderaxespad'] = -1.6
plt.rcParams['axes.linewidth']    = 1
mpl.rcParams['lines.linewidth']   = 2

# Define color palette and colormap
palette = {'heatmap' : 'viridis',
           'black'   : '#000000',
           'blue'    : '#377eb8',
           'green'   : '#4daf4a',
           'red'     : '#e41a1c'}

def plot_window(title, ylabel, xlabel):
    """Initializes a new plotting window."""

    plt.figure(tight_layout = True)
    plt.gcf().canvas.set_window_title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

def plot_density(data):
    """Function for plotting the density of states."""

    # Extract and reshape data
    position = np.unique(data[:,0])
    energy   = np.unique(data[:,1])
    density  = data[:,2].reshape(len(position), len(energy))

    # Plot the density of states
    plot_window('Density of states', r'Position $z/\xi$', r'Energy $\epsilon/\Delta$')
    plt.pcolormesh(energy, position, density, vmin = 0.0, vmax = 2.0, cmap = palette['heatmap'])
    plt.colorbar(ticks = [0.0,0.5,1.0,1.5,2.0])
    plt.axis([-2.0, +2.0, position.min(), position.max()])
    plt.show()

def plot_current(data):
    """Function for plotting the charge and spin currents."""

    # Extract and reshape data
    position = np.unique(data[:,0])
    charge   = data[:,1]
    spin_x   = data[:,2]
    spin_y   = data[:,3]
    spin_z   = data[:,4]

    # Plot the charge current
    plot_window('Charge current', r'Charge current $I_{e}/I_{e0}$', r'Position $z/\xi$')
    plt.plot(position, charge, palette['black'])
    plt.axis([position.min(), 
              position.max(), 
              min(0.0, charge.min()) * 1.05,
              max(0.0, charge.max()) * 1.05])
    plt.show()

    # Plot the spin current
    plot_window('Spin current', r'Spin current $I_{\sigma}/I_{\sigma 0}$', r'Position $z/\xi$')
    plt.plot(position, spin_x, palette['blue'], 
             position, spin_y, palette['green'], 
             position, spin_z, palette['red'])
    plt.axis([position.min(), 
              position.max(), 
              min(0.0, spin_x.min(), spin_y.min(), spin_z.min())*1.05, 
              max(0.0, spin_x.max(), spin_y.max(), spin_z.max())*1.05])
    plt.legend([r'Spin-$x$', r'Spin-$y$', r'Spin-$z$'], ncol = 3)
    plt.show()

def plot_gap(data):
    """Function for plotting the superconducting gap and phase."""

    # Extract and reshape data
    position = np.unique(data[:,0])
    gap      = data[:,1]
    phase    = data[:,2]

    # Plot the superconducting gap
    plt.figure(tight_layout = True)
    plt.gcf().canvas.set_window_title('Superconducting gap')
    plt.plot(position, gap, palette['black'])
    plt.axis([position.min(), 
              position.max(), 
              min(-1e-6, gap.min()) * 1.05,
              max(+1e-6, gap.max()) * 1.05])
    plt.ylabel(r'Superconducting gap $\Delta/\Delta_0$')
    plt.xlabel(r'Position $z/\xi$')
    plt.show()

    # Plot the superconducting phase
    plt.figure(tight_layout = True)
    plt.gcf().canvas.set_window_title('Superconducting phase')
    plt.plot(position, phase, palette['black'])
    plt.axis([position.min(), 
              position.max(), 
              min(-1e-6, phase.min()) * 1.05,
              max(+1e-6, phase.max()) * 1.05])
    plt.ylabel(r'Superconducting phase $\varphi/\pi$')
    plt.xlabel(r'Position $z/\xi$')
    plt.show()

if __name__ == "__main__":
    # Check command line arguments
    try:
        filename = sys.argv[1]
    except:
        print('Usage: {} <filename.dat>'.format(sys.argv[0]))
        sys.exit(0)

    # Read data from file
    with open(filename, 'r') as f:
        desc = f.readline()   # File description
        data = np.loadtxt(f)  # File contents

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
