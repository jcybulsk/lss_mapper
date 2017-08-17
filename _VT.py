import numpy as np
import pandas as pd
import matplotlib
import scipy
import seaborn as sns

import pylab as pl
from scipy.spatial import Voronoi
import matplotlib.pyplot as plt
from matplotlib import rc
import imageio
import warnings; warnings.filterwarnings('ignore', message='elementwise')
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import (FixedLocator, MaxNLocator,
                                                 DictFormatter)
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes

#### MST clustering library and astropy for coordinate transforms
sns.set_style("whitegrid")

""" Compute area of a polynomial given input vertices """
def get_polyarea(x, y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

class VT:
    """ Voronoi Tessellation Class 
    
    Compute Voronoi tessellation on a set of input points using 
    scipy's Voronoi function.
    
    ---- Input ----
    data frame "df", which has two columns present:
     - ra: right ascension (deg)
     - dec: declination (deg)
     
    ---- Attributes ----
    celldensity: base-10 log of surface density for each input 
    galaxy (in units of 1/[sq. deg]), calculated from the inverse 
    of the cell area
    
    """
    def __init__(self, df, celldensity=None):
        self.df = df
        self.celldensity = celldensity
        self.get_vt(df)

    """ Run scipy Voronoi function """
    def get_vt(self, df):
        pos = np.array([list(i) for i in zip(df.ra, df.dec)])
        return Voronoi(pos)
        
    """ Calculate cell densities """
    def get_celldens(self, V):
        holdx = []
        holdy = []
        celldens = []
        for reg in V.regions:
            if -1 in reg:
                celldens.append(-99)
            else:
                for num in reg:
                    holdx.append(V.vertices[num][0])
                    holdy.append(V.vertices[num][1])
                cellarea = get_polyarea(holdx, holdy)
                if cellarea == 0:
                    # There will always be one entry in the 
                    # Voronoi regions list of lists what will be 
                    # and empty list []. That's this case. It 
                    # should be left out, so just do nothing here.
                    pass
                else:
                    celldens.append(np.log10(1.0/float(cellarea)))
                holdx = []
                holdy = []
        return celldens
    
    """ Visualize the Voronoi cells """
    def plot_voronoi(self, V, *args, **kwargs):
        xlim = kwargs.get('xlim', None)
        ylim = kwargs.get('ylim', None)
        savefigure = kwargs.get('savefigure', False)
        figname = kwargs.get('figname', 'VT_figure.png')
        plt.figure(figsize=(14,8))
        for vind in V.ridge_vertices:
            (i1, i2) = sorted(vind)
            if (i1 != -1) & (i2 != -1):
                vor1 = V.vertices[i1]
                vor2 = V.vertices[i2]
                p, = plt.plot([vor1[0], vor2[0]], [vor1[1], vor2[1]], 'k-')
        plt.xlabel('Right Ascension (deg)', fontsize=18)
        plt.ylabel('Declination (deg)', fontsize=18)
        if xlim != None:
            plt.xlim(xlim)
        if ylim != None:
            plt.ylim(ylim)

        # Leave an option to save all the plots to output PNG files.
        if savefigure == True:
        	pl.savefig(figname, bbox_inches='tight', dpi=250)
            
    """ Plot the cumulative distribution of VT cell densities """
    def plot_vt_cumul(self, V, *args, **kwargs):
        savefigure = kwargs.get('savefigure', False)
        figname = kwargs.get('figname', 'VT_cumul_dist.png')
        vcalc = np.asarray(V.celldensity)
        celldens_calc = vcalc[vcalc > -99]
        sns.distplot(celldens_calc, hist_kws=dict(cumulative=False), 
                     kde_kws=dict(cumulative=True))
        plt.xlabel('log$_{10}$ (VT Cell Densities)', fontsize=15)
        plt.ylabel('Norm. Counts/Cumul. Dist.', fontsize=15)
        
        # Leave an option to save all the plots to output PNG files.
        if savefigure == True:
        	pl.savefig(figname, bbox_inches='tight')