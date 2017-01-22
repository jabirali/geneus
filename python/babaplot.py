import string
import matplotlib as mpl
import matplotlib.pyplot as plt

# Set the matplotlib stylesheet
mpl.style.use('~/Code/python/babaplot.mplstyle')

# Define functions for creating figures
def onecolumn(n, m):
    '''Return a figure that is one column wide.'''
    return plt.subplots(n, m, figsize=(1*3.35, 1*(2.2*m)/n+0.3*(m>1)), tight_layout=True)

def twocolumn(n, m):
    '''Return a figure that are two columns wide.'''
    return plt.subplots(n, m, figsize=(2*3.35, 2*(2.2*m)/n+0.3*(m>1)), tight_layout=True)

# Define functions for styling figures
def compactify(fig, label=1):
    '''Restyle a figure to save space. This is appropriate for one-column figures with multiple plots stacked next to each other horizontally.'''
    for n, ax in enumerate(fig.axes):
        # Label each subfigure in order
        if (label > 0):
            # Calculate the label placement
            if (label == 1):
                va = 'top'
                ha = 'right'
                xy = (0.98, 0.96)
            elif (label == 2):
                va = 'top'
                ha = 'left'
                xy = (0.02, 0.96)
            elif (label == 3):
                va = 'bottom'
                ha = 'left'
                xy = (0.02, 0.04)
            else:
                va = 'bottom'
                ha = 'right'
                xy = (0.98, 0.04)

            # Generate the label
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
