import numpy as np
import pandas as pd
import matplotlib
import seaborn as sns

import pylab as pl
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
from mst_clustering import MSTClustering
from astropy.coordinates import SkyCoord
import astropy.units as u
sns.set_style("whitegrid")

class MST:
    """ Minimum Spanning Tree Class
    
    Compute MST for a set of input points using the MSTClustering 
    code from jakevdp, calculate branch lengths from the MST and 
    generate plots of the MST and cumulative distribution of branch 
    lengths.
    
    ---- Inputs ----
    data frame "df", which has two columns present:
     - ra: right ascension (deg)
     - dec: declination (deg)
    
    cutoff_scale (float): minimum size of edges, also known as the 
                          critical branch length. All edges larger 
                          than cutoff_scale will be removed.
    
    min_cluster_size (int): min number of galaxies in a cluster.
    
    n_neighbors (int): maximum number of neighbors of each point 
    used for approximate Euclidean MST algorithm.
    
    ---- Attributes ----
    labels: integer specifying the structure to which a given galaxy 
            has been assigned. It will have a -1 if no membership was 
            assigned.
            
    segments: sets of ra, dec coordinates for the MST branch segments
    seps: base-10 log of branch lengths (in degrees)
    
    """
    def __init__(self, df, cutoff_scale=None, min_cluster_size=None, 
                 n_neighbors=None, set_mst=None, labels=None, 
                 segments=None, seps=None):
        self.df = df
        self.cutoff_scale = cutoff_scale
        self.min_cluster_size = min_cluster_size
        self.n_neighbors = n_neighbors
        self.set_mst = MSTClustering(cutoff_scale=cutoff_scale, 
                             min_cluster_size=min_cluster_size, 
                             n_neighbors=n_neighbors)
        pos = np.array([list(i) for i in zip(df.ra, df.dec)])
        self.labels = self.set_mst.fit_predict(pos)
        self.segments = self.set_mst.get_graph_segments(full_graph=True)
        self.seps = self.get_sep_mst()
        
    """ Calculate branch lengths (in base-10 log(degrees)) 
        from the MST segments """    
    def get_sep_mst(self):
        mst_coord0_ra = np.asarray(self.segments[0][0])
        mst_coord1_ra = np.asarray(self.segments[0][1])
        mst_coord0_dec = np.asarray(self.segments[1][0])
        mst_coord1_dec = np.asarray(self.segments[1][1])
        c0 = SkyCoord(mst_coord0_ra, mst_coord0_dec, unit=u.deg)
        c1 = SkyCoord(mst_coord1_ra, mst_coord1_dec, unit=u.deg)
        return np.log10(c0.separation(c1).degree)
        
    """ Plot the MST diagram (left) and the labeled structures 
        identified from the MST (right) """
    def plot_mst(self, model, cmap='rainbow', *args, **kwargs):
        """Utility code to visualize a minimum spanning tree"""
        xlim = kwargs.get('xlim', None)
        ylim = kwargs.get('ylim', None)
        ssize = kwargs.get('s', 8)
        savefigure = kwargs.get('savefigure', False)
        X = model.X_fit_
        fig, ax = plt.subplots(1, 2, figsize=(20, 7), sharex=True, sharey=True)
        for axi, full_graph, colors in zip(ax, [True, False], ['lightblue', model.labels_]):
            segments = model.get_graph_segments(full_graph=full_graph)
            axi.plot(segments[0], segments[1], '-k', zorder=1, lw=1)
            plt.xlabel('Right Ascension (deg)', size=14)
            plt.ylabel('Declination (deg)', size=14)
            axi.scatter(X[:, 0], X[:, 1], c=colors, cmap=cmap, zorder=2, s=ssize)
            axi.axis('tight')
            if xlim != None:
                plt.xlim(xlim)
            if ylim != None:
                plt.ylim(ylim)
    
        ax[0].set_title('Full Minimum Spanning Tree', size=16)
        ax[1].set_title('Trimmed Minimum Spanning Tree', size=16)
        
        # Leave an option to save all the plots to output PNG files.
        if savefigure == True:
        	pl.savefig('MST_figure.png', bbox_inches='tight')
        
        
    """ Plot the cumulative distribution of MST branch lengths """
    def plot_mst_cumul(self, *args, **kwargs):
    	savefigure = kwargs.get('savefigure', False)
        sns.distplot(self.seps, hist_kws=dict(cumulative=False), 
                     kde_kws=dict(cumulative=True))
        plt.xlabel('log$_{10}$ (MST branch length)', fontsize=15)
        plt.ylabel('Norm. Counts/Cumul. Dist.', fontsize=15)

        # Leave an option to save all the plots to output PNG files.
        if savefigure == True:
        	pl.savefig('MST_cumul_dist.png', bbox_inches='tight')