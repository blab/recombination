'''
This script produces a heatmap showing genetic distance between each segment for random pairs of samples.
The heatmap is ordered by full genome genetic distance.

Inputs are:
    --pairwise, pickle file containing comparison dictionary
    --max-distance, maximum y of plot.
    --output, name of output PNG
'''

import argparse
import pickle
import random
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def load(pfile):
    '''
    Loads pairwise dictionary from pickle file.
    '''
    with open(pfile, 'rb') as file:
        pairwise = pickle.load(file)
    return pairwise

def choose_random_pairs(pairwise, max_distance):
    '''
    Chooses a random pair of samples with a full-genome genetic distance of X for every X in range(0,max_distance)
    '''
    iterate_list = [*range(max_distance)]
    pair_list = []
    for num in iterate_list:
        distance = []
        for pair in pairwise:
            if pairwise[pair]['distance'] == num:
                distance.append(pair)
        pair_list.append(distance)
    pairs = []
    for distance in pair_list:
        pairs.append(random.choice(distance))
    return pairs

def pairs_to_array(max_distance, segments, pair_list, pairwise):
    '''
    Creates a matrix where each row represents is a pairwise, and each column is the genetic distance between each segment.
    Rows are ordered from 0 to max_distance.
    '''
    array = np.zeros((max_distance, len(segments)), dtype=int)
    for pair, row in zip(pair_list, range(max_distance)):
        for segment, value in zip(segments, range(len(segments))):
            array[row, value] = pairwise[pair][segment]['distance']
    return array

def heatmap(array, segments, max_distance, output):
    '''
    Plots array as a heatmap.
    '''
    mpl.rcParams['font.weight']=200
    mpl.rcParams['axes.labelweight']=200
    mpl.rcParams['font.size']=14

    fig, ax = plt.subplots(figsize=(8, 10), facecolor='white')
    image = ax.pcolormesh(array)
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label('Genetic distance\nper segment')

    ax.set_xticks(np.arange(len(segments))+0.5)
    ax.set_yticklabels(range(0, max_distance, 50))
    ax.set_xticklabels(segments, ha = 'center')
    ax.set_ylabel('Pairwise genetic distance')

    return plt.savefig(output, dpi=250)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Heatmap of genetic distance for influenza',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--pairwise', type=str, required=True, help='pairwise pickle file')
    parser.add_argument('--max-distance', type=int, required=True, help='Maximum distance for heatmap')
    parser.add_argument('--output', type=str, required=True, help = 'name of output figure')
    args = parser.parse_args()

    # Define influenza segments
    segments = ['ha', 'na', 'pb2', 'pb1', 'pa', 'np', 'mp', 'ns']

    # Loads comparison dictionary from pickle file
    pairwise = load(args.pairwise)

    # Chooses random pairs for every distance from 0, max_distance
    random_pairs = choose_random_pairs(pairwise, args.max_distance)

    # Creates array to plot.
    array = pairs_to_array(args.max_distance, segments, random_pairs, pairwise)

    # Plots heatmap of array
    heatmap(array, segments, args.max_distance, args.output)
