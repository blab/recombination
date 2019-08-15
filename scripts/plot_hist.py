'''
This script plots histograms of pairwise divergence for the full influenza genome and each influenza genome segment.
Histograms are colored for within cluster pairs and between cluster pairs for clusters based on specific genetic distance cutoff.

Inputs are:
    --pairwise, a pickle dictionary with gentic distance for sample pairs.
    --cutoff, Hamming distance at which to cluster samples
    --output-genome, name of PNG for full genome histogram
    --output-segments, name of PNG for segments histogram
'''

import argparse
import pickle
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import random
from scipy import stats
import matplotlib as mpl
from matplotlib import pyplot as plt

def strains(pairwise):
    '''
    Returns set of the strains compared in pairwise dictionary.
    '''
    strain_list = set()
    for (strainA, strainB) in pairwise.keys():
        strain_list.add(strainA)
        strain_list.add(strainB)
    return strain_list

def adjacency_matrix(strain_list, pairwise, segments, cutoff):
    '''
    Returns an adjacency matrix on which to cluster.
    '''
    adj_matrix = np.zeros((len(strain_list), len(strain_list)))
    for indexA, strainA in enumerate(strain_list):
        for indexB, strainB in enumerate(strain_list):
            if indexA != indexB:
                if (strainA, strainB) in pairwise:
                    distance = pairwise[(strainA, strainB)]['distance']
                elif (strainB, strainA) in pairwise:
                    distance = pairwise[(strainB, strainA)]['distance']
                if distance < cutoff:
                    adj_matrix[indexA, indexB] = 1
                    adj_matrix[indexB, indexA] = 1
    return adj_matrix

def cluster(adj_matrix, strain_list):
    '''
    Uses SciPy connected_components to cluster samples.
    '''
    graph = csr_matrix(adj_matrix)
    n_clusters, labels = connected_components(
        csgraph=graph, directed=False, return_labels=True
    )
    cluster_dict =  dict(zip(strain_list, labels))
    return cluster_dict

def zscore(data):
    '''
    Returns z-score.
    '''
    z = np.abs(stats.zscore(data))
    return z

def make_distances_dict(pairwise, segments, cluster_dict):
    '''
    Creates dictionary containing key : list of genetic distances for 10,000 random pairs of samples.
    Keys include 'all', 'between_cluster', 'within_cluster', and within & between for each influenza segment.
    '''
    keys = random.sample(pairwise.keys(), 10000)
    distances = {}
    distances['all'] = []
    distances['within'] = []
    distances['between'] = []
    for segment in segments:
        distances[segment + '_within'] = []
        distances[segment + '_between'] = []
    for key in keys:
        distance = pairwise[key]['distance']
        distances['all'].append(distance)
        (strainA, strainB) = key
        if cluster_dict[strainA] == cluster_dict[strainB]:
            distances['within'].append(distance)
        else:
            distances['between'].append(distance)
        for segment in segments:
            segment_distance = pairwise[key][segment]['distance']
            if cluster_dict[strainA] == cluster_dict[strainB]:
                distances[segment + '_within'].append(segment_distance)
            else:
                distances[segment + '_between'].append(segment_distance)
    return distances

def remove_outliers(distance_dict):
    '''
    Removes outliers from distances to plot.
    '''
    distances_edited = {}
    for key, value in distance_dict.items():
        z = zscore(value)
        distances_edited[key] = [distance for index, distance in enumerate(distance_dict[key]) if z[index] < 4]
    return distances_edited

def hist_genome(distance_dict, output_genome):
    '''
    Plots histogram of genetic distance for entire flu genome
    '''
    mpl.rcParams['font.weight']=110
    mpl.rcParams['axes.labelweight']=110
    mpl.rcParams['font.size']=14

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 12), facecolor='white', sharex=True, sharey=True)
    ax1.set_title('H3N2: Full genome')
    ax1.hist(distance_dict['all'], bins=len(set(distance_dict['all'])))
    ax1.set_ylabel('Frequency')
    ax2.hist((distance_dict['within'], distance_dict['between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set(distance_dict['all'])))
    ax2.set_xlabel('Pairwise genetic distance')
    ax2.set_ylabel('Frequency')
    ax2.legend()
    return fig.savefig(output_genome, dpi=300)

def hist_segments(distance_dict, output_segments):
    '''
    Plots histogram of genetic distance for each individual flu segment.
    '''
    mpl.rcParams['font.weight']=110
    mpl.rcParams['axes.labelweight']=110
    mpl.rcParams['font.size']=14

    fig, axs = plt.subplots(2, 4, figsize=(16, 12), facecolor='white', sharex=True, sharey=True)
    fig.suptitle('H3N2', fontsize='large')
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
    ax2.set_title('NA')
    ax2.hist((distance_dict['na_within'], distance_dict['na_between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set.union(set(distance_dict['na_within']), set(distance_dict['na_between']))))
    ax3.set_title('PB2')
    ax3.hist((distance_dict['pb2_within'], distance_dict['pb2_between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set.union(set(distance_dict['pb2_within']), set(distance_dict['pb2_between']))))
    ax4.set_title('PB1')
    ax4.hist((distance_dict['pb1_within'], distance_dict['pb1_between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set.union(set(distance_dict['pb1_within']), set(distance_dict['pb1_between']))))
    ax5.set_title('PA')
    ax5.hist((distance_dict['pa_within'], distance_dict['pa_between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set.union(set(distance_dict['pa_within']), set(distance_dict['pa_between']))))
    ax6.set_title('NP')
    ax6.hist((distance_dict['np_within'], distance_dict['np_between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set.union(set(distance_dict['np_within']), set(distance_dict['np_between']))))
    ax7.set_title('MP')
    ax7.hist((distance_dict['mp_within'], distance_dict['mp_between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set.union(set(distance_dict['mp_within']), set(distance_dict['mp_between']))))
    ax8.set_title('NS')
    ax8.hist((distance_dict['ns_within'], distance_dict['ns_between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set.union(set(distance_dict['ns_within']), set(distance_dict['ns_between']))))

    return fig.savefig(output_segments, dpi=300)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plot distribution of genetic distance.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--pairwise', type=str, required=True, help='pairwise pickle file')
    parser.add_argument('--cutoff', type=int, required=True, help='Genomic distance cutoff on which to cluster')
    parser.add_argument('--output-genome', type=str, required=True, help = 'name of output figure')
    parser.add_argument('--output-segments', type=str, required=True, help = 'name of output figure')
    args = parser.parse_args()

    # Define influenza segments
    segments = ['ha', 'na', 'pb2', 'pb1', 'pa', 'np', 'mp', 'ns']

    # Loads dictionar(y/ies) containing pairwise divegence
    with open(args.pairwise, 'rb') as file:
        pairwise = pickle.load(file)

    # Lists strains shared across all segments
    strains_list = strains(pairwise)

    # Defines adjacency matrix on which to cluster
    adj_matrix = adjacency_matrix(strains_list, pairwise, segments, args.cutoff)

    # Clusters strains based on genetic distance
    clusters = cluster(adj_matrix, strains_list)

    # Creates list of distances to plot
    distances = make_distances_dict(pairwise, segments, clusters)

    # Removes outliers with a Z-score > 4 from distances dictionary.
    distances_edited = remove_outliers(distances)

    # Writes histograms of distances
    hist_genome(distances_edited, args.output_genome)
    hist_segments(distances_edited, args.output_segments)
