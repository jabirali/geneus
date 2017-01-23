import string
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# Set the pandas output format
pd.options.display.width = 132

# Set the matplotlib stylesheet
mpl.style.use('~/Code/python/babaplot.mplstyle')



################################################################################
# Subroutines for loading and processing data
################################################################################

def load(filename, statements=[], index=0):
    '''Define a function for loading and interpolating a dataframe.'''

    # Read the file into a dataframe
    df = pd.read_table(filename, delim_whitespace=True, index_col=index)

    # Remove any duplicate entries
    df.drop_duplicates(inplace=True)

    # Process any eval statements
    for statement in statements:
        df.eval(statement, inplace=True)

    # Generate a new index
    index = np.linspace(df.index.min(), df.index.max(), 1024)

    # Absorb the old index
    index = np.unique(np.concatenate((df.index, index)))

    # Reindex and interpolate
    return df.reindex(index).interpolate(method='pchip')



################################################################################
# Subroutines for plotting data
################################################################################

def grid(cols, n, m):
    '''Return a figure that is one column wide.'''
    return plt.subplots(n, m, figsize=(cols*3.35, cols*(2.2*n)/m+0.3*(m>1)), tight_layout=True)

def compactify(fig):
    '''Restyle a figure to save space. This is appropriate for one-column figures with multiple plots stacked next to each other horizontally.'''
    for n, ax in enumerate(fig.axes):
        # Label each subfigure in order
        va = 'bottom'
        ha = 'center'
        xy = (0.5, 1.04)

        ax.annotate('(' + string.ascii_lowercase[n] + ')', xy=xy, va=va, ha=ha, 
                xycoords='axes fraction', fontweight='bold', fontsize='x-small')

        # Disable the legend
        ax.get_legend().set_visible(False)

        # Reposition the axis labels
        ax.get_yaxis().set_label_coords(-0.05, 0.50)
        ax.get_xaxis().set_label_coords( 0.50,-0.12)

        # Only display ticks at the endpoints
        ax.set_xticks(ax.get_xlim())
        ax.set_yticks(ax.get_ylim())

        # Use scientific notation on the yticks
        ax.ticklabel_format(style='sci', axis='x', scilimits=(-3,+3))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(-0,+0))
