#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import itertools
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Define plot style
plt.style.use('~/Code/python/babaplot.mplstyle')

def plot_window(xlabel, ylabel, title):
    """Initializes a new plotting window."""

    plt.figure(tight_layout = True)
    plt.gcf().canvas.set_window_title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    return plt.gca()

def multiplot(lst, title='Plot'):
    '''Function for plotting multiple plots in the same window, specified as a list of funtions.'''

    # Create a figure window
    fig = plt.figure(tight_layout = True)
    fig.canvas.set_window_title(title)

    def onclick(event):
        fig.clear()
        next(onclick.lst)()
        fig.canvas.draw()

    onclick.lst = itertools.cycle(lst)
    onclick(None)
    fig.canvas.mpl_connect('button_press_event', onclick)

    plt.show()

def plot_density(data):
    """Function for plotting the density of states."""

    # Extract and reshape data
    position = np.unique(data[:,0])
    energy   = np.unique(data[:,1])
    density  = data[:,2].reshape(len(position), len(energy))

    # Plot the density of states
    ax = plot_window(r'Energy $ε/Δ$', r'Position $z/ξ$', 'Density of states')
    plt.pcolormesh(energy, position, density, vmin = 0.0, vmax = 2.0)
    plt.colorbar(ticks = [0.0,0.5,1.0,1.5,2.0])
    plt.axis([-2.0, +2.0, position.min(), position.max()])
    plt.show()

def plot_magnetization(data):
    """Function for plotting the magnetization."""

    # Extract and reshape data
    position = np.unique(data[:,0])
    mag_x   = data[:,1]
    mag_y   = data[:,2]
    mag_z   = data[:,3]

    # Plot the data
    ax = plot_window(r'Position $z/ξ$', r'Magnetization $M/M_0$', 'Magnetization')
    plt.plot(position, mag_x, 
             position, mag_y, 
             position, mag_z)
    plt.axis([position.min(), 
              position.max(), 
              min(0.0, mag_x.min(), mag_y.min(), mag_z.min())*1.05, 
              max(0.0, mag_x.max(), mag_y.max(), mag_z.max())*1.05])
    plt.legend([r'$M_x$', r'$M_y$', r'$M_z$'], ncol = 3)
    plt.show()

def plot_current(data):
    """Function for plotting the charge and spin currents."""

    # Extract and reshape data
    position = np.unique(data[:,0])
    charge   = data[:,1]
    spin_x   = data[:,2]
    spin_y   = data[:,3]
    spin_z   = data[:,4]

    # Plot the spin current
    ax = plot_window(r'Position $z/ξ$', r'Spin current $I_{σ}/I_{σ0}$', 'Spin current')
    plt.plot(position, spin_x, 
             position, spin_y, 
             position, spin_z)
    plt.axis([position.min(), 
              position.max(), 
              min(0.0, spin_x.min(), spin_y.min(), spin_z.min())*1.05, 
              max(0.0, spin_x.max(), spin_y.max(), spin_z.max())*1.05])
    plt.legend([r'Spin-$x$', r'Spin-$y$', r'Spin-$z$'], ncol = 3)

    # Plot the charge current
    ax = plot_window(r'Position $z/ξ$', r'Charge current $I_{e}/I_{e0}$', 'Charge current')
    plt.plot(position, charge)
    plt.axis([position.min(), 
              position.max(), 
              min(0.0, charge.min()) * 1.05,
              max(0.0, charge.max()) * 1.05])

    plt.show()

def plot_gap(data):
    """Function for plotting the superconducting gap and phase."""

    # Extract and reshape data
    position = np.unique(data[:,0])
    gap      = data[:,1]
    phase    = data[:,2]

    def plot_magnitude():
        # Plot the superconducting gap
        plt.plot(position, gap)
        plt.axis([position.min(), 
                  position.max(), 
                  min(-1e-6, gap.min()) * 1.05,
                  max(+1e-6, gap.max()) * 1.05])
        plt.ylabel(r'Superconducting gap $\Delta/\Delta_0$')
        plt.xlabel(r'Position $z/\xi$')

    def plot_phase():
        # Plot the superconducting phase
        plt.plot(position, phase)
        plt.axis([position.min(), 
                  position.max(), 
                  min(-1e-6, phase.min()) * 1.05,
                  max(+1e-6, phase.max()) * 1.05])
        plt.ylabel(r'Superconducting phase $\varphi/\pi$')
        plt.xlabel(r'Position $z/\xi$')

    multiplot([plot_magnitude, plot_phase], title = 'Superconducting gap')


if __name__ == "__main__":
    # Check command line arguments
    if len(sys.argv) == 1:
        print('Usage: {} <filename.dat>'.format(sys.argv[0]))
        sys.exit(0)

    # Loop over input files
    for filename in sys.argv[1:]:
        print(filename)

        # Read data from file
        with open(filename, 'r') as f:
            desc = f.readline()   # File description
            data = np.loadtxt(f)  # File contents

        if filename.endswith('.dat') or filename.endswith('.bak'):
            # Split the file description into fields
            field = [q.strip() for q in desc[2:].split('\t')]

            # Call the appropriate plotting method
            if field == ['Position', 'Energy', 'Density of states']:
                plot_density(data)
            elif field == ['Position', 'Charge current', 'Spin-x current', 'Spin-y current', 'Spin-z current']:
                plot_current(data)
            elif field == ['Position', 'Gap magnitude', 'Gap phase']:
                plot_gap(data)
            #elif field == ['Position', 'Mx magnetization', 'My magnetization', 'Mz magnetization']:
            else:
                plot_magnetization(data)
            #else:
            #    print('Sorry, unrecognized output file...')
        elif filename.endswith('.conf'):
            print('Sorry, conf parser not integrated yet...')
        else:
            print('Sorry, unrecognized data format...')
