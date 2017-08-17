import numpy as np
import matplotlib.pyplot as plt
import pylab as pl

def plot_both_cumul_dist(cell_dens, mst_branches, *args, **kwargs):
    savefigure = kwargs.get('savefigure', False)
    dens_thresh = kwargs.get('dens_thresh', None)
    cell_dens = sorted(cell_dens[cell_dens > -99])
    mst_branches = sorted(mst_branches)
    cell_dens_cumul = 1.0*np.arange(len(cell_dens))/(len(cell_dens) - 1)
    mst_branches_cumul = 1.0*np.arange(len(mst_branches))/(len(mst_branches) - 1)

    x_min = min([min(cell_dens), min(mst_branches)])
    x_max = max([max(cell_dens), max(mst_branches)])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(cell_dens, cell_dens_cumul, color='b')
    plt.plot(mst_branches, mst_branches_cumul, color='r')
    plt.xlim(x_min, x_max)
    plt.ylim(0, 1)
    
    for i in range(len(dens_thresh)):
        dens_interp = np.interp(dens_thresh[i], cell_dens, cell_dens_cumul)
        plt.axvline(x=dens_thresh[i], ymin=0, ymax=dens_interp, 
                    color='blue', linestyle='--')
        plt.hlines(y=dens_interp, xmin=x_min, xmax=dens_thresh[i], 
                    color='blue', linestyle='--')
        
        mst_cumul_val = 1.0-dens_interp
        mst_interp = np.interp(mst_cumul_val, mst_branches_cumul, mst_branches)
        plt.vlines(x=mst_interp, ymin=0, ymax=mst_cumul_val, 
                    color='red', linestyle='--')
        plt.hlines(y=mst_cumul_val, xmin=x_min, xmax=mst_interp, 
                    color='red', linestyle='--')
                    
    ax.text(x_min+0.3, 0.8, 'VT', fontsize=15, color='blue')
    ax.text(x_min+0.3, 0.7, 'MST', fontsize=15, color='red')
    plt.xlabel('log$_{10}$ (VT Cell Dens or MST branch len)', fontsize=15)
    plt.ylabel('Cumulative Fraction', fontsize=14)
    
    # Leave an option to save all the plots to output PNG files.
    if savefigure == True:
    	pl.savefig('both_cumul_dist.png', bbox_inches='tight')