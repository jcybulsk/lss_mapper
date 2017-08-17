import numpy as np
import pandas as pd
import matplotlib
import seaborn as sns

import pylab as pl
import matplotlib.pyplot as plt
from matplotlib import rc
import imageio

sns.set_style("whitegrid")

def plot_both_MST(model, labels1, labels2, *args, **kwargs):
    xlim = kwargs.get('xlim', None)
    ylim = kwargs.get('ylim', None)
    ssize = kwargs.get('s', 8)
    savefigure = kwargs.get('savefigure', False)
    figname = kwargs.get('figname', 'MST_combined.png')
    cmap = "rainbow"
        
    # To get the colors of the points right, first set the filament structures 
    # to a specific index integer (like 50?), and then set all the 
    labels_overall = labels1
    labels_overall[labels1 > -1] = 50
    labels_overall[labels2 > -1] = 150
    X = model.X_fit_
    
    colors = labels_overall
    segments = model.get_graph_segments(full_graph=False)
    plt.figure(figsize=(10,6))
    plt.plot(segments[0], segments[1], '-k', zorder=1, lw=1)
    plt.xlabel('Right Ascension (deg)', size=14)
    plt.ylabel('Declination (deg)', size=14)
    plt.scatter(X[:, 0], X[:, 1], c=colors, cmap=cmap, zorder=2, s=ssize)
    plt.axis('tight')
    if xlim != None:
        plt.xlim(xlim)
    if ylim != None:
        plt.ylim(ylim)
        
    # Leave an option to save all the plots to output PNG files.
    if savefigure == True:
    	pl.savefig(figname, bbox_inches='tight', dpi=250)
