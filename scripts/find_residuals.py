'''
This script finds the cumulative residuals for pairwise TMRCA between all influenza
segment combinations. It stores those residuals in a NxN matrix in the HDF5 file
and produces a histogram of their distribution.
'''

import argparse
import h5py
import numpy as np
from itertools import combinations
import random
import matplotlib as mpl
import matplotlib.pyplot as plt

def transform_to_tmrca(segments, file, clock_rate, segment_length):
    '''
    Returns dictionary containing NxN matrices of TMRCA for each influenza segment.
    '''
    tmrca = {}
    for segment, rate, length in zip(segments, clock_rate, segment_length):
        matrix = np.array(file['samples'][segment].get(segment))
        segment_rate = rate * length
        matrix = (matrix / segment_rate)
        tmrca[segment] = matrix
    return tmrca

def calculate_residuals(segments, tmrca):
    '''
    Returns NxN matrix containing cumulative residuals for pairwise TMRCA between all influenza segments.
    '''
    residuals = np.zeros((tmrca['ha'].shape))
    combos = combinations(segments, 2)
    for (segmentA, segmentB) in combos:
        residuals += np.abs(tmrca[segmentA] - tmrca[segmentB])
    return residuals

def choose_random_indices(residuals):
    '''
    Chooses 10k random pairs to plot.
    '''
    indices = np.tril_indices(residuals.shape[0], -1)
    keys = random.sample(list(zip(*indices)), 10000)
    keys_0 = np.asarray([i for i,j in keys])
    keys_1 = np.asarray([j for i, j in keys])
    return keys_0, keys_1

def scatter_tmrca(lineage, tmrca, keys_0, keys_1, output):
    '''
    Plots TMRCA vs. TMRCA of various segments for 10k random pairs.
    '''
    mpl.rcParams['font.weight']=110
    mpl.rcParams['axes.labelweight']=110
    mpl.rcParams['font.size']=14
    mpl.rcParams["scatter.marker"]='.'
    mpl.rcParams['lines.markersize']=2
    mpl.rcParams['lines.linewidth'] = 0.7

    fig, (axs) = plt.subplots(2, 2, figsize=(10, 8), facecolor='white', sharex=True, sharey=True)
    fig.suptitle(lineage, fontsize='large')
    ax1 = axs[0,0]
    ax2 = axs[0,1]
    ax3 = axs[1,0]
    ax4 = axs[1,1]

    ax1.scatter(x = tmrca['ha'][keys_0, keys_1], y = tmrca['na'][keys_0, keys_1])
    ax1.set_xlabel('TMRCA (years) in HA')
    ax1.set_ylabel('TMRCA (years) in NA')
    ax1.set_xlim(left=0, right=13)
    ax1.set_ylim(bottom=0, top=13)
    x = np.linspace(*ax1.get_xlim())
    ax1.plot(x,x, c='k')

    ax2.scatter(x = tmrca['ha'][keys_0, keys_1], y = tmrca['pb2'][keys_0, keys_1])
    ax2.set_xlabel('TMRCA (years) in HA')
    ax2.set_ylabel('TMRCA (years) in PB2')
    ax2.plot(x,x, c='k')

    ax3.scatter(x = tmrca['na'][keys_0, keys_1], y = tmrca['pb2'][keys_0, keys_1])
    ax3.set_xlabel('TMRCA (years) in NA')
    ax3.set_ylabel('TMRCA (years) in PB2')
    ax3.plot(x,x, c='k')

    ax4.scatter(x = tmrca['pa'][keys_0, keys_1], y = tmrca['pb1'][keys_0, keys_1])
    ax4.set_xlabel('TMRCA (years) in PA')
    ax4.set_ylabel('TMRCA (years) in PB1')
    ax4.plot(x,x, c='k')
    return fig.savefig(output, dpi=300)

def plot_residuals(residuals, keys_0, keys_1, lineage, output):
    '''
    Plots histogram of residuals for the 10k random pairs.
    '''
    residuals_sample = residuals[keys_0, keys_1]

    mpl.rcParams['font.weight']=110
    mpl.rcParams['axes.labelweight']=110
    mpl.rcParams['font.size']=14
    bins = [n for n in range(0, 250, 5)]

    fig, ax = plt.subplots(figsize=(10, 8), facecolor='white')
    ax.set_title(lineage)
    ax.set_xlim(left=0, right=250)
    ax.hist(residuals_sample, bins=bins, align='mid', edgecolor='k')
    ax.set_xlabel('TMRCA (years) across all segments')
    ax.set_ylabel('Frequency')
    return fig.savefig(output, dpi=300)

def write_to_h5py(file, matrix):
    '''
    Writes residuals matrix to HDF5 file under ['samples']['residuals']
    '''
    group = file['samples']['genome']
    if 'residuals' in group:
        del group['residuals']
    group.create_dataset('residuals', data=matrix)
    return file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Find TMRCA residuals',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--pairwise', type=str, required=True, help='pairwise hdf5 file')
    parser.add_argument('--clock-rate', type=float, nargs='+', required=True, help='list of clock rates for lineage')
    parser.add_argument('--segment-length', type=int, nargs='+', required=True, help='list of segment lengths for lineage')
    parser.add_argument('--lineage', type=str, help = 'name of viral lineage')
    parser.add_argument('--output-scatter', type=str, required=True, help = 'name of output scatter plots')
    parser.add_argument('--output-hist', type=str, required=True, help = 'name of output histogram')
    args = parser.parse_args()

    # Define influenza segments
    segments = ['ha', 'na', 'pb2', 'pb1', 'pa', 'np', 'mp', 'ns']

    # Opens HDF5 file
    hfile = h5py.File(args.pairwise, mode='r+')

    # Transforms genetic distance into TMRCA
    tmrca = transform_to_tmrca(segments, hfile, args.clock_rate, args.segment_length)

    # Calculates cumulative residuals
    residuals = calculate_residuals(segments, tmrca)

    # Chooses 10K random pairs
    keys0, keys1 = choose_random_indices(residuals)

    # Plots TMRCA vs. TMRACA for various segments
    scatter_tmrca(args.lineage, tmrca, keys0, keys1, args.output_scatter)

    # Plots histogram of cumulative residuals
    plot_residuals(residuals, keys0, keys1, args.lineage, args.output_hist)

    # Saves to HDF5 file and closes HDF5 file
    write_to_h5py(hfile, residuals)
