'''
This scripts finds pairwise divergence for each sample pair in an alignment.
Currently, it works for influenza alignments given in the order: 'ha', 'na', 'pb2', 'pb1' 'pa', 'np', 'mp', 'ns'.

Its inputs are:
    --alignments, a list of genome segment alignments
    --output, the location of an output dictionary.
'''

import argparse
from Bio import SeqIO
import numpy as np
import h5py

def find_shared(align_list):
    '''
    Returns list of strains shared across all genome segments.
    '''
    setlist = []
    for alignment in align_list:
        strains = set()
        for record in SeqIO.parse(alignment, 'fasta'):
            strains.add(record.id)
        setlist.append(strains)
    shared = set.intersection(*setlist)
    return shared

def convert(align_list, segments, shared):
    '''
    Replaces nucleotide bases with numbers in order to find SNPs efficiently. C,T,G,A are replaced with the first 4 prime numbers.
    Ambiguous nucleotide bases are replaced with the multiple of those primes.
    For example, R, which represents A or G, is 35 because A is 7, and G is 5.
    '''
    mapping = {}
    for strain in shared:
        mapping[strain] = {}
    for alignment, segment in zip(align_list, segments):
        for record in SeqIO.parse(alignment, 'fasta'):
            if record.id in mapping:
                seq = list(record.seq)
                seq_conversion = {
                    'C': 2,
                    'T': 3,
                    'G': 5,
                    'A': 7,
                    'D': 105,
                    'R': 35,
                    'Y': 6,
                    'H': 42,
                    'K': 15,
                    'V': 70,
                    'N': 210,
                    'M': 14,
                    'W': 21,
                    'S': 10
                }
                sequence = [seq_conversion.get(base) for base in seq]
                array = np.asarray(sequence, dtype=int)
                mapping[record.id][segment] = array
    return mapping

def compare(array1, array2):
    '''
    Identifies SNPs between sequences where sequences are represented as numpy integer arrays.
    CTGA must be represented by prime numbers and ambiguous bases are the multiples of those prime numbers.
    '''
    condition1 = (array1 != array2)
    condition2 = (array1 % array2 != 0)
    condition3 = (array2 % array1 != 0)
    location = np.where(condition1 & condition2 & condition3)
    distance = np.size(location)
    return distance

def compare_pairwise(shared, mapping, segments):
    '''
    Compares all strains pairwise storing the genetic distance between each strains as NxN arrays at the genome & segment level.
    '''
    counter = 0
    interval = 500
    length = len(shared)
    pairwise = {}
    for segment in segments:
        matrix = np.zeros((length, length), dtype='i4')
        for indexA, strainA in enumerate(shared):
            if counter % interval == 0:
                print("[", end = '')
                for x in range(int(counter/interval)):
                    print("-", end = '')
                for x in range(int((length*len(segments))/interval) - int(counter/interval)):
                    print(" ", end = '')
                print("]")
            for indexB, strainB in enumerate(shared):
                matrix[indexA, indexB] = compare(mapping[strainA][segment], mapping[strainB][segment])
            counter +=1
        pairwise[segment] = np.tril(matrix, -1)
    pairwise['genome'] = np.zeros((length, length))
    for segment in segments:
        pairwise['genome'] = np.add(pairwise['genome'], pairwise[segment])
    return pairwise

def write_to_h5py(output, shared, pairwise):
    with h5py.File(output, mode='a') as file:
        samples = file.create_group('samples')
        shared_strains = np.asarray(shared, dtype='a50')
        samples.create_dataset('samples', data=shared_strains)
        for key in pairwise:
            subgrp = samples.create_group(key)
            subgrp.create_dataset(key, data=pairwise[key], compression='gzip')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Pairwise analysis for each sequence.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--alignments', type=str, nargs='+', required=True, help='alignment for viruses')
    parser.add_argument('--output', type=str, required=True, help = 'name of output dictionary')
    args = parser.parse_args()

    # Define influenza segments.
    segments = ['ha', 'na', 'pb2', 'pb1', 'pa', 'np', 'mp', 'ns']

    # Identifies set of strains shared across all segments.
    shared_strains = find_shared(args.alignments)

    # Converts bases from letters to numbers for comparison
    converted = convert(args.alignments, segments, shared_strains)

    # Calculates pairwise divergence for every strain
    pairwise = compare_pairwise(shared_strains, converted, segments)

    # Saves matrices to HDF5 file.
    write_to_h5py(args.output, shared_strains, pairwise)
