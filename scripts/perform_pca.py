'''
This script uses PCA to visualize viral strain clusters to assess if PCA
can be used to cluster influenza constellations:

(1) Performs PCA on sample pairs, coloring by within and between cluster pairs.

(2) Performs PCA on samples, coloring by clusters. Matrix rows are samples in N.
    Matrix columns are genetic distance between each segment for every strain in N.
'''

import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sklearn
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def create_matrix_1(file, segments, cutoff):
    '''
    Returns a 2D array where each row is a pairwise comparison. Columns are the hamming
    distance between each segment, ordered according to 'ha', 'na', 'pb2', 'pb1','pa', 'np', 'mp', 'ns'.
    Also returns a 1D array ordered, according to matrix row pairs. If the pairs are in the same cluster, the array = 1. Otherwise, array = 0.
    '''
    array = np.array(file['samples']['genome'].get('genome'))
    indices = np.tril_indices(array.shape[0], -1)
    indices_list = list(zip(*indices))
    matrix = np.zeros((len(indices_list), len(segments)))
    for segment, column in zip(segments, range(len(segments))):
        segment_matrix = np.array(file['samples'][segment].get(segment))
        for index, row in zip(indices_list, range(len(indices_list))):
            matrix[row, column] = segment_matrix[index]

    clusters = np.array(file['samples'].get('clusters'+ str(cutoff)))
    strains_0 = np.asarray([i for i,j in indices_list])
    strains_1 = np.asarray([j for i, j in indices_list])
    clusters_0 = clusters[strains_0]
    clusters_1 = clusters[strains_1]
    cluster_status = np.zeros((matrix.shape[0],))
    cluster_status[clusters_0 == clusters_1] = 1
    return matrix, cluster_status

def create_matrix_2(file, segments, cutoff):
    '''
    Returns N x (N*8) array. Each row is a sample. Columns are genetic distance between influenza segments for all samples.
    '''
    array = np.array(file['samples']['genome'].get('genome'))
    matrix = np.zeros((array.shape[0], (array.shape[0]*8)))
    counter = 0
    for segment in segments:
        segment_arr = np.array(file['samples'][segment].get(segment))
        arr_full = np.triu(segment_arr.T) + segment_arr
        n_col = arr_full.shape[1]
        matrix[:, counter : n_col+counter] = arr_full
        counter += n_col

    clusters = np.array(file['samples'].get('clusters' + str(cutoff)))
    return matrix, clusters

def standardize_data(array):
    '''
    Standardizes segment distances for PCA.
    '''
    array_standardized = StandardScaler().fit_transform(array)
    return array_standardized

def perform_pca(array, components):
    '''
    Returns principal components matrix.
    '''
    pca = PCA(n_components = components, svd_solver='full')
    principal_components = pca.fit_transform(array)
    return principal_components, pca.explained_variance_ratio_

def add_target(array, target):
    '''
    Adds targets for coloring to the principal components matrix.
    '''
    array_final = np.zeros((array.shape[0], array.shape[1] + 1))
    array_final[0:array.shape[0], 0:array.shape[1]] = array
    array_final[:, -1] = target.T
    return array_final

def plot_pca_1(principal_components, lineage, explained_variance, output):
    '''
    Plots explained variance & PC1, PC2, & PC3 for PCA (1)
    '''
    array = principal_components[np.random.choice(principal_components.shape[0], 10000, replace=False), :]

    mpl.rcParams['font.weight']=110
    mpl.rcParams['axes.labelweight']=110
    mpl.rcParams['font.size']=14
    mpl.rcParams["scatter.marker"]='.'
    mpl.rcParams['lines.markersize']=2

    fig, (axs) = plt.subplots(2, 2, figsize=(20, 16), facecolor='white')
    fig.suptitle(lineage, fontsize='large')
    ax1 = axs[0, 0]
    ax2 = axs[0, 1]
    ax3 = axs[1, 0]
    ax4 = axs[1, 1]

    ax1.plot(explained_variance)
    ax1.set_xlabel('Principal components')
    ax1.set_ylabel('Ratio of explained variance')

    for cluster in np.unique(array[:,-1]):
        if cluster == 1:
            label = 'Within cluster'
        elif cluster == 0:
            label = 'Between cluster'
        pc1 = array[array[:,-1] == cluster, 0]
        pc2 = array[array[:,-1] == cluster, 1]
        pc3 = array[array[:,-1] == cluster, 2]

        ax2.scatter(x=pc1, y=pc2, label=label)
        ax2.set_xlabel('Principal Component 1')
        ax2.set_ylabel('Principal Component 2')
        ax2.legend()

        ax3.scatter(x=pc1, y=pc3, label=label)
        ax3.set_xlabel('Principal Component 1')
        ax3.set_ylabel('Principal Component 3')
        ax3.legend()

        ax4.scatter(x=pc2, y=pc3, label=label)
        ax4.set_xlabel('Principal Component 2')
        ax4.set_ylabel('Principal Component 3')
        ax4.legend()
    return plt.savefig(output, dpi=300)

def plot_pca_2(array, lineage, explained_variance, output):
    '''
    Plots explained variance & PC1, PC2, & PC3 for PCA (2)
    '''
    mpl.rcParams['font.weight']=110
    mpl.rcParams['axes.labelweight']=110
    mpl.rcParams['font.size']=14
    mpl.rcParams["scatter.marker"]='.'
    mpl.rcParams['lines.markersize']=6

    fig, (axs) = plt.subplots(2, 2, figsize=(20, 16), facecolor='white')
    fig.suptitle(lineage, fontsize='large')
    ax1 = axs[0, 0]
    ax2 = axs[0, 1]
    ax3 = axs[1, 0]
    ax4 = axs[1, 1]

    ax1.plot(explained_variance)
    ax1.set_xlabel('Principal components')
    ax1.set_ylabel('Ratio of explained variance')

    for cluster in np.unique(array[:,-1]):
        pc1 = array[array[:,-1] == cluster, 0]
        pc2 = array[array[:,-1] == cluster, 1]
        pc3 = array[array[:,-1] == cluster, 2]

        ax2.scatter(x=pc1, y=pc2, label=cluster)
        ax2.set_xlabel('Principal Component 1')
        ax2.set_ylabel('Principal Component 2')

        ax3.scatter(x=pc1, y=pc3, label=cluster)
        ax3.set_xlabel('Principal Component 1')
        ax3.set_ylabel('Principal Component 3')

        ax4.scatter(x=pc3, y=pc2, label=cluster)
        ax4.set_xlabel('Principal Component 3')
        ax4.set_ylabel('Principal Component 2')
    return plt.savefig(output, dpi=300)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Perform PCA on strains and pairs of strains',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--pairwise', type=str, required=True, help='pairwise hdf5 file')
    parser.add_argument('--cutoff', type=int, required=True, help='Cutoff for clustering')
    parser.add_argument('--components', type=int, default=10, help='No. of components for PCA')
    parser.add_argument('--lineage', type=str, required=True, help = 'name of viral lineage')
    parser.add_argument('--output-pairs', type=str, required=True, help = 'name of figure for PCA (1)')
    parser.add_argument('--output-strains', type=str, required=True, help = 'name of figure for PCA (2)')
    args = parser.parse_args()

    # Define influenza segments
    segments = ['ha', 'na', 'pb2', 'pb1', 'pa', 'np', 'mp', 'ns']

    # Opens HDF5 file
    hfile = h5py.File(args.pairwise, mode='r')

    # Creates matrix for PCA (1)
    matrix_1, target_1 = create_matrix_1(hfile, segments, args.cutoff)

    # Creates matrix for PCA (2)
    matrix_2, target_2 = create_matrix_2(hfile, segments, args.cutoff)

    # Standardizes matrices for PCA
    standardized_1 = standardize_data(matrix_1)
    standardized_2 = standardize_data(matrix_2)

    # Performs PCA
    pc_matrix_1, exp_var_1 = perform_pca(standardized_1, args.components)
    pc_matrix_2, exp_var_2 = perform_pca(standardized_2, args.components)

    # Adds targets to principal component matrices
    matrix_to_plot_1 = add_target(pc_matrix_1, target_1)
    matrix_to_plot_2 = add_target(pc_matrix_2, target_2)

    # Plot PCA (1)
    plot_pca_1(matrix_to_plot_1, args.lineage, exp_var_1, args.output_pairs)

    # Plot PCA (2)
    plot_pca_2(matrix_to_plot_2, args.lineage, exp_var_2, args.output_strains)

    # Closes HDF5 File
    hfile.close()
