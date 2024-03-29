'''
This script produces a heatmap showing genetic distance between each segment for random pairs of samples.
The heatmap is ordered by full genome genetic distance.

Inputs are:
    --pairwise, HDF5 database produced by compare.py
    --max-distance, maximum y of plot.
    --clock-rate, list specifying molecular clock rate for each influenza segment
    --segment-length, list specifying the length of each influenza segment
    --lineage, lineage of virus, e.g. 'h3n2'
    --output, name of output PNG
'''

import argparse
import h5py
import random
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def choose_random_pairs(file, max_distance):
    '''
    Chooses a random pair of samples with a full-genome genetic distance of X for every X in range(0,max_distance)
    '''
    matrix = np.array(file['samples']['genome'].get('genome'))
    lower_tri = np.tril(np.ones(matrix.shape), -1)
    pair_list = []
    for i in range(max_distance):
        pair_list.append(np.where((lower_tri == 1) & (matrix == i)))
    random_pairs = []
    for pairs in pair_list:
        random_pairs.append(random.choice(list(zip(*pairs))))
    return random_pairs

def pairs_to_array(max_distance, segments, pair_list, file, clock_rate, segment_length):
    '''
    Creates a matrix where each row represents a pair of strains, and each column is the genetic distance between each segment.
    Rows are ordered from 0 to max_distance.
    '''
    array = np.zeros((max_distance, len(segments)), dtype=int)
    for pair, row in zip(pair_list, range(max_distance)):
        for segment, value in zip(segments, range(len(segments))):
            array[row, value] = np.array(file['samples'][segment].get(segment))[pair]
    segment_rate = np.asarray([rate * length for rate, length in zip(clock_rate, segment_length)])
    array = (array / segment_rate)
    return array

def heatmap(array, segments, max_distance, lineage, output):
    '''
    Plots array as a heatmap.
    '''
    mpl.rcParams['font.weight']=200
    mpl.rcParams['axes.labelweight']=200
    mpl.rcParams['font.size']=14

    fig, ax = plt.subplots(figsize=(8, 10), facecolor='white')
    image = ax.pcolormesh(array, vmax=25)
    cbar = fig.colorbar(image, ax=ax, extend='max')
    cbar.set_label('Years since divergence')

    ax.set_xticks(np.arange(len(segments))+0.5)
    ax.set_yticklabels(range(0, max_distance, 50))
    ax.set_xticklabels(segments, ha = 'center')
    ax.set_ylabel('Pairwise genetic distance')
    ax.set_title(lineage)
    return plt.savefig(output, dpi=300)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Heatmap of genetic distance for influenza',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--pairwise', type=str, required=True, help='pairwise hdf5 file')
    parser.add_argument('--max-distance', type=int, required=True, help='Maximum distance for heatmap')
    parser.add_argument('--clock-rate', type=float, nargs='+', required=True, help='list of clock rates for lineage')
    parser.add_argument('--segment-length', type=int, nargs='+', required=True, help='list of segment lengths for lineage')
    parser.add_argument('--lineage', type=str, help = 'name of viral lineage')
    parser.add_argument('--output', type=str, required=True, help = 'name of output figure')
    args = parser.parse_args()

    # Define influenza segments
    segments = ['ha', 'na', 'pb2', 'pb1', 'pa', 'np', 'mp', 'ns']

    # Opens HDF5 file
    hfile = h5py.File(args.pairwise, mode='r')

    # Chooses random pairs for every distance from 0, max_distance
    random_pairs = choose_random_pairs(hfile, args.max_distance)

    # Creates array to plot.
    array = pairs_to_array(args.max_distance, segments, random_pairs, hfile, args.clock_rate, args.segment_length)

    # Plots heatmap of array
    heatmap(array, segments, args.max_distance, args.lineage, args.output)

    # Closes HDF5 file
    hfile.close()
