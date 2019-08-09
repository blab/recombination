
import argparse
import pickle
import numpy as np
from itertools import combinations
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import random
from scipy import stats
import matplotlib as mpl
from matplotlib import pyplot as plt

def cluster(shared, pairwise_dict, cutoff):
    adj_matrix = np.zeros((len(shared), len(shared)))
    for indexA, strainA in enumerate(shared):
        for indexB, strainB in enumerate(shared):
            if indexA != indexB:
                if (strainA, strainB) in pairwise_dict:
                    if pairwise_dict[(strainA, strainB)]['distance'] < cutoff:
                        adj_matrix[indexA, indexB] = 1
                        adj_matrix[indexB, indexA] = 1
                elif pairwise_dict[(strainB, strainA)]['distance'] < cutoff:
                    adj_matrix[indexA, indexB] = 1
                    adj_matrix[indexB, indexA] = 1

    graph = csr_matrix(adj_matrix)
    n_clusters, labels = connected_components(
        csgraph=graph, directed=False, return_labels=True
    )

    cluster_dict =  dict(zip(shared, labels))
    return cluster_dict

def zscore(data):
    z = np.abs(stats.zscore(data))
    return z

def list_distances(pairwise_dict, cluster_dict):
    keys = list(random.sample(pairwise_dict.keys(), 10000))
    distances = {}
    distances['all_dist'] = []
    distances['within_cluster'] = []
    distances['between_cluster'] = []
    distances['ha'] = []
    distances['ha_within'] = []
    distances['ha_between'] = []
    distances['na'] = []
    distances['na_within'] = []
    distances['na_between'] = []
    distances['pb2'] = []
    distances['pb2_within'] = []
    distances['pb2_between'] = []
    distances['pb1'] = []
    distances['pb1_within'] = []
    distances['pb1_between'] = []
    distances['pa'] = []
    distances['pa_within'] = []
    distances['pa_between'] = []
    distances['np'] = []
    distances['np_within'] = []
    distances['np_between'] = []
    distances['mp'] = []
    distances['mp_within'] = []
    distances['mp_between'] = []
    distances['ns'] = []
    distances['ns_within'] = []
    distances['ns_between'] = []

    for key in keys:
        distances['all_dist'].append(pairwise_dict[key]['distance'])
        ha_location = [loc for loc in pairwise_dict[key]['snp_loc'] if loc < 1702]
        ha_distance = np.size(ha_location)
        distances['ha'].append(ha_distance)
        na_location = [loc for loc in pairwise_dict[key]['snp_loc'] if loc >= 1702 & loc < 3138]
        na_distance = np.size(na_location)
        distances['na'].append(na_distance)
        pb2_location = [loc for loc in pairwise_dict[key]['snp_loc'] if loc >= 3138 & loc < 5448]
        pb2_distance = np.size(pb2_location)
        distances['pb2'].append(pb2_distance)
        pb1_location = [loc for loc in pairwise_dict[key]['snp_loc'] if loc >=5448 & loc < 7759]
        pb1_distance = np.size(pb1_location)
        distances['pb1'].append(pb1_distance)
        pa_location = [loc for loc in pairwise_dict[key]['snp_loc'] if loc >= 7759 & loc < 9951]
        pa_distance = np.size(pa_location)
        distances['pa'].append(pa_distance)
        np_location = [loc for loc in pairwise_dict[key]['snp_loc'] if loc >= 9951 & loc < 11488]
        np_distance = np.size(np_location)
        distances['np'].append(np_distance)
        mp_location = [loc for loc in pairwise_dict[key]['snp_loc'] if loc >= 11488 & loc < 12487]
        mp_distance = np.size(mp_location)
        distances['mp'].append(mp_distance)
        ns_location = [loc for loc in pairwise_dict[key]['snp_loc'] if loc >= 12487]
        ns_distance = np.size(ns_location)
        distances['ns'].append(ns_distance)
        (strainA, strainB) = key
        if cluster_dict[strainA] == cluster_dict[strainB]:
            distances['within_cluster'].append(pairwise_dict[key]['distance'])
            distances['ha_within'].append(ha_distance)
            distances['na_within'].append(na_distance)
            distances['pb2_within'].append(pb2_distance)
            distances['pb1_within'].append(pb1_distance)
            distances['pa_within'].append(pa_distance)
            distances['np_within'].append(np_distance)
            distances['mp_within'].append(mp_distance)
            distances['ns_within'].append(ns_distance)
        else:
            distances['between_cluster'].append(pairwise_dict[key]['distance'])
            distances['ha_between'].append(ha_distance)
            distances['na_between'].append(na_distance)
            distances['pb2_between'].append(pb2_distance)
            distances['pb1_between'].append(pb1_distance)
            distances['pa_between'].append(pa_distance)
            distances['np_between'].append(np_distance)
            distances['mp_between'].append(mp_distance)
            distances['ns_between'].append(ns_distance)

    # Removs outliers from data
    distances_edited = {}
    for key, value in distances.items():
        z = zscore(value)
        distances_edited[key] = [distance for index, distance in enumerate(distances[key]) if z[index] < 4]

    return distances_edited


def list_distances_ha(pairwise_dict, cluster_dict):
    keys = list(random.sample(pairwise_dict.keys(), 10000))
    distances = {}
    distances['all_dist'] = []
    distances['within_cluster'] = []
    distances['between_cluster'] = []
    for key in keys:
        distances['all_dist'].append(pairwise_dict[key]['distance'])
        (strainA, strainB) = key
        if cluster_dict[strainA] == cluster_dict[strainB]:
            distances['within_cluster'].append(pairwise_dict[key]['distance'])
        else:
            distances['between_cluster'].append(pairwise_dict[key]['distance'])

    # Removes outliers from data. There ended up being one outlier at a genetic distance of 1750, that I think might just be the wrong sequence.
    distances_edited = {}
    z_all = zscore(distances['all_dist'])
    z_within = zscore(distances['within_cluster'])
    z_between = zscore(distances['between_cluster'])
    distances_edited['all_dist'] = [distance for index, distance in enumerate(distances['all_dist']) if z_all[index] < 4]
    distances_edited['within_cluster'] = [distance for index, distance in enumerate(distances['within_cluster']) if z_within[index] < 4]
    distances_edited['between_cluster'] = [distance for index, distance in enumerate(distances['between_cluster']) if z_between[index] < 4]

    return distances_edited

def hist_distances(distance_dict, output_genome, output_segments):
    mpl.rcParams['font.weight']=110
    mpl.rcParams['axes.labelweight']=110
    mpl.rcParams['font.size']=14

    fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 12), facecolor="white")
    ax1.set_title('H3N2: Full genome')
    ax1.hist(distance_dict['all_dist'], bins=len(set(distance_dict['all_dist'])))
    ax1.set_ylabel("Frequency")
    ax2.hist((distance_dict['within_cluster'], distance_dict['between_cluster']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set(distance_dict['all_dist'])))
    ax2.set_xlabel("Pairwise genetic distance")
    ax2.set_ylabel("Frequency")
    ax2.legend()

    fig2, axs = plt.subplots(4, 4, figsize=(16, 12), facecolor="white", sharex=True, sharey=True)
    fig2.suptitle('H3N2', fontsize='large')
    fig2.text(0.5, 0.04, 'Pairwise genetic distance', ha='center')
    fig2.text(0.04, 0.5, 'Frequency', va='center', rotation='vertical')
    ax1 = axs[0,0]
    ax2 = axs[1,0]
    ax3 = axs[0,1]
    ax4 = axs[1,1]
    ax5 = axs[0,2]
    ax6 = axs[1,2]
    ax7 = axs[0,3]
    ax8 = axs[1,3]
    ax9 = axs[2,0]
    ax10 = axs[3,0]
    ax11 = axs[2,1]
    ax12 = axs[3,1]
    ax13 = axs[2,2]
    ax14 = axs[3,2]
    ax15 = axs[2,3]
    ax16 = axs[3,3]

    ax1.hist(distance_dict['ha'], label='HA', bins=len(set(distance_dict['ha'])))
    ax1.set_title('HA')
    ax2.hist((distance_dict['ha_within'], distance_dict['ha_between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set(distance_dict['ha'])))
    ax3.hist(distance_dict['na'], label='NA', bins=len(set(distance_dict['na'])))
    ax3.set_title('NA')
    ax4.hist((distance_dict['na_within'], distance_dict['na_between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set(distance_dict['na'])))
    ax5.hist(distance_dict['pb2'], label='PB2', bins=len(set(distance_dict['pb2'])))
    ax5.set_title('PB2')
    ax6.hist((distance_dict['pb2_within'], distance_dict['pb2_between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set(distance_dict['pb2'])))
    ax7.hist(distance_dict['pb1'], label='PB1', bins=len(set(distance_dict['pb1'])))
    ax7.set_title('PB1')
    ax8.hist((distance_dict['pb1_within'], distance_dict['pb1_between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set(distance_dict['pb1'])))
    ax9.hist(distance_dict['pa'], label='PA', bins=len(set(distance_dict['pa'])))
    ax9.set_title('PA')
    ax10.hist((distance_dict['pa_within'], distance_dict['pa_between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set(distance_dict['pa'])))
    ax11.hist(distance_dict['np'], label='NP', bins=len(set(distance_dict['np'])))
    ax11.set_title('NP')
    ax12.hist((distance_dict['np_within'], distance_dict['np_between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set(distance_dict['np'])))
    ax13.hist(distance_dict['mp'], label='MP', bins=len(set(distance_dict['mp'])))
    ax13.set_title('MP')
    ax14.hist((distance_dict['mp_within'], distance_dict['mp_between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set(distance_dict['mp'])))
    ax15.hist(distance_dict['ns'], label='NS', bins=len(set(distance_dict['ns'])))
    ax15.set_title('NS')
    ax16.hist((distance_dict['ns_within'], distance_dict['ns_between']), label = ('Within clusters', 'Between clusters'), stacked=True, bins=len(set(distance_dict['ns'])))

    return fig1.savefig(output_genome, dpi=300), fig2.savefig(output_segments, dpi=300)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Plot distribution of genetic distance.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--shared', type=str, required=True, help="name of shared strains pickle file")
    parser.add_argument('--pairwise', type=str, required=True, help="name of pairwise pickle file")
    parser.add_argument('--cutoff', type=int, required=True, help="Genomic distance cutoff on which to cluster")
    parser.add_argument('--output-genome', type=str, required=True, help = "name of output figure")
    parser.add_argument('--output-segments', type=str, required=True, help = "name of output figure")
    args = parser.parse_args()

    # Loads list of shared strains
    with open(args.shared, "rb") as file:
        shared = pickle.load(file)

    # Loads dictionary containing pairwise divegence
    with open(args.pairwise, "rb") as file:
        pairwise = pickle.load(file)

    # Cluster strains by genetic distance
    clusters = cluster(shared, pairwise, args.cutoff)

    # Creates list of distances to plot
    distances = list_distances(pairwise, clusters)

    # Writes histogram of distances
    hist_distances(distances, args.output_genome, args.output_segments)
