'''
This script plots histograms of pairwise divergence for the full influenza genome and each influenza genome segment.
Histograms are colored for within cluster pairs and between cluster pairs for clusters based on a specific genetic distance cutoff.

Inputs are:
    --pairwise, HDF5 file created by compare.py.
    --cutoff, Hamming distance at which to cluster samples
    --lineage, lineage of virus, e.g. 'h3n2'
    --output-genome, name of PNG for full genome histogram
    --output-segments, name of PNG for segments histogram
'''

import argparse
import h5py
import numpy as np
import random
from scipy import stats
import matplotlib as mpl
from matplotlib import pyplot as plt

def make_distances_dict(file, segments, cutoff):
    '''
    Creates dictionary containing key : list of genetic distances for 10,000 random pairs of samples.
    Keys include 'all', 'between_cluster', 'within_cluster', and within & between for each influenza segment.
    '''
    matrix = np.array(file['samples']['genome'].get('genome'))
    indices = np.tril_indices(matrix.shape[0], -1)
    keys = random.sample(list(zip(*indices)), 10000)
    clusters = np.array(file['samples'].get('clusters' + str(cutoff)))

    keys_0 = np.asarray([i for i,j in keys])
    keys_1 = np.asarray([j for i, j in keys])
    clusters_0 = clusters[keys_0]
    clusters_1 = clusters[keys_1]
    loc_within = np.where(clusters_0 == clusters_1)
    loc_between = np.where(clusters_0 != clusters_1)
    within_keys_0 = keys_0[loc_within]
    within_keys_1 = keys_1[loc_within]
    between_keys_0 = keys_0[loc_between]
    between_keys_1 = keys_1[loc_between]

    distances = {}
    distances['all'] = matrix[keys_0, keys_1]
    distances['within'] = matrix[within_keys_0, within_keys_1]
    distances['between'] = matrix[between_keys_0, between_keys_1]

    for segment in segments:
        segment_matrix = np.array(file['samples'][segment].get(segment))
        distances[segment + '_within'] = segment_matrix[within_keys_0, within_keys_1]
        distances[segment + '_between'] = segment_matrix[between_keys_0, between_keys_1]
    return distances

def zscore(data):
    '''
    Returns z-score.
    '''
    z = np.abs(stats.zscore(data))
    return z

def remove_outliers(distance_dict):
    '''
    Removes outliers from distances to plot.
    '''
    distances_edited = {}
    for key, value in distance_dict.items():
        z = zscore(value)
        distances_edited[key] = [distance for index, distance in enumerate(distance_dict[key]) if z[index] < 4]
    return distances_edited

def hist_genome(distance_dict, lineage, output_genome):
    '''
    Plots histogram of genetic distance for entire flu genome
    '''
    mpl.rcParams['font.weight']=110
    mpl.rcParams['axes.labelweight']=110
    mpl.rcParams['font.size']=14

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 16), facecolor='white', sharex=True, sharey=True)
    ax1.set_title(lineage+ ': full genome')
    ax1.hist(distance_dict['all'], bins=len(set(distance_dict['all'])))
    ax1.set_ylabel('Frequency')
    ax1.set_xlim(left=0)
    ax2.hist((distance_dict['within'], distance_dict['between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set(distance_dict['all'])))
    ax2.set_xlabel('Pairwise genetic distance')
    ax2.set_ylabel('Frequency')
    ax2.set_xlim(left=0)
    ax2.legend()
    return fig.savefig(output_genome, dpi=300)

def hist_segments(distance_dict, lineage, output_segments):
    '''
    Plots histogram of genetic distance for each individual flu segment.
    '''
    mpl.rcParams['font.weight']=110
    mpl.rcParams['axes.labelweight']=110
    mpl.rcParams['font.size']=14

    fig, axs = plt.subplots(2, 4, figsize=(12, 16), facecolor='white', sharex=True, sharey=True)
    fig.suptitle(lineage, fontsize='large')
    fig.text(0.5, 0.04, 'Pairwise genetic distance', ha='center')
    fig.text(0.04, 0.5, 'Frequency', va='center', rotation='vertical')
    ax1 = axs[0,0]
    ax2 = axs[0,1]
    ax3 = axs[0,2]
    ax4 = axs[0,3]
    ax5 = axs[1,0]
    ax6 = axs[1,1]
    ax7 = axs[1,2]
    ax8 = axs[1,3]

    ax1.set_title('HA')
    ax1.hist((distance_dict['ha_within'], distance_dict['ha_between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set.union(set(distance_dict['ha_within']), set(distance_dict['ha_between']))))
    ax1.set_xlim(left=0)
    ax2.set_title('NA')
    ax2.hist((distance_dict['na_within'], distance_dict['na_between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set.union(set(distance_dict['na_within']), set(distance_dict['na_between']))))
    ax2.set_xlim(left=0)
    ax3.set_title('PB2')
    ax3.hist((distance_dict['pb2_within'], distance_dict['pb2_between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set.union(set(distance_dict['pb2_within']), set(distance_dict['pb2_between']))))
    ax3.set_xlim(left=0)
    ax4.set_title('PB1')
    ax4.hist((distance_dict['pb1_within'], distance_dict['pb1_between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set.union(set(distance_dict['pb1_within']), set(distance_dict['pb1_between']))))
    ax4.set_xlim(left=0)
    ax5.set_title('PA')
    ax5.hist((distance_dict['pa_within'], distance_dict['pa_between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set.union(set(distance_dict['pa_within']), set(distance_dict['pa_between']))))
    ax5.set_xlim(left=0)
    ax6.set_title('NP')
    ax6.hist((distance_dict['np_within'], distance_dict['np_between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set.union(set(distance_dict['np_within']), set(distance_dict['np_between']))))
    ax6.set_xlim(left=0)
    ax7.set_title('MP')
    ax7.hist((distance_dict['mp_within'], distance_dict['mp_between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set.union(set(distance_dict['mp_within']), set(distance_dict['mp_between']))))
    ax7.set_xlim(left=0)
    ax8.set_title('NS')
    ax8.hist((distance_dict['ns_within'], distance_dict['ns_between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set.union(set(distance_dict['ns_within']), set(distance_dict['ns_between']))))
    ax8.set_xlim(left=0)
    return fig.savefig(output_segments, dpi=300)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plot distribution of genetic distance.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--pairwise', type=str, required=True, help='HDF5 dataset generated by compare.py')
    parser.add_argument('--cutoff', type=int, required=True, help='Genomic distance cutoff on which to cluster')
    parser.add_argument('--lineage', type=str, required=True, help='lineage of virus')
    parser.add_argument('--output-genome', type=str, required=True, help = 'name of output figure')
    parser.add_argument('--output-segments', type=str, required=True, help = 'name of output figure')
    args = parser.parse_args()

    # Define influenza segments
    segments = ['ha', 'na', 'pb2', 'pb1', 'pa', 'np', 'mp', 'ns']

    # Opens HDF5 file
    hfile = h5py.File(args.pairwise, mode='r')

    # Creates list of distances to plot
    distances = make_distances_dict(hfile, segments, args.cutoff)

    # Removes outliers with a Z-score > 4 from distances dictionary.
    distances_edited = remove_outliers(distances)

    # Writes histograms of distances
    hist_genome(distances_edited, args.lineage, args.output_genome)
    hist_segments(distances_edited, args.lineage, args.output_segments)

    # Closes HDF5 file
    hfile.close()
