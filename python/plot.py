#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import configparser as cp
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# Set the matplotlib stylesheet
mpl.style.use('babaplot.mplstyle')

# Read and parse a config file
config = cp.ConfigParser()
config.read('plot.conf')

files = {filename: columns.replace(' ','').split(',') for filename, columns in config['files'].items()}
plots = config['plots']
axes  = config['axes']

# Read and concatenate input files into one dataframe
df = pd.concat((pd.read_table(key, names=files[key]) for key in sorted(files)), axis=1)

# Plot the datasets
for label, data in plots.items():
    # Define the expression to plot
    x_data, y_data = data.replace(' ','').split(',')
    df.eval('_x = ' + x_data, inplace=True)
    df.eval('_y = ' + y_data, inplace=True)

    # Perform the actual plotting
    plt.plot('_x', '_y', data=df, label=label)

# Finalize
plt.xlabel(axes.get('x', ''))
plt.ylabel(axes.get('y', ''))
plt.legend(ncol=5)
plt.show()
