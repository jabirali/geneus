#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import configparser
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.style.use('babaplot.mplstyle')


# What to read and plot -- get this from user somehow
ax_lab = { 'x': 'Horizontal axis',
           'y': 'Vertical axis'  }

# Read config file
config = configparser.ConfigParser()
config.read('plot.conf')

files = {filename: columns.replace(' ','').split(',') for filename, columns in config['files'].items()}

x_data = config['plots']['x'].replace(' ','')
y_data = config['plots']['y'].replace(' ','').split(',')

# Read and concatenate input files into one dataframe
df = pd.concat((pd.read_table(key, names=files[key]) for key in sorted(files)), axis=1)

# Define the expressions that will be plotted
df.eval('_x = ' + x_data, inplace=True)

# Perform the plotting using matplotlib
for y in y_data:
    df.eval('_y = ' + y, inplace=True)
    plt.plot('_x', '_y', data=df)
    #df.eval('_y = ' + y, inplace=True)
    #plt.plot('_x', '_y', data=df, label=label)

# Finalize
plt.xlabel(ax_lab['x'])
plt.ylabel(ax_lab['y'])
plt.legend(ncol=len(y_data))
plt.show()
