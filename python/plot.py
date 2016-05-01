#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.style.use('babaplot.mplstyle')

# What to read and plot -- get this from user somehow
files = {'critical_0.000.dat': ['a', 'b', 'c'],
         'critical_0.500.dat': ['u', 'v', 'w']}

ax_lab = { 'x': 'Horizontal axis',
           'y': 'Vertical axis'  }

x_data = 'a'
y_data = { 'c':   'c data',
           'c-w': 'c-w data'}

# Read and concatenate input files into one dataframe
df = pd.concat((pd.read_table(key, names=files[key]) for key in sorted(files)), axis=1)

# Define the expressions that will be plotted
df.eval('_x = ' + x_data, inplace=True)

# Perform the plotting using matplotlib
for y, label in y_data.items():
    df.eval('_y = ' + y, inplace=True)
    plt.plot('_x', '_y', data=df, label=label)

# Finalize
plt.xlabel(ax_lab['x'])
plt.ylabel(ax_lab['y'])
plt.legend(ncol=len(y_data))
plt.show()
