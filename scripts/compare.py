'''
This script takes an alignment of sequences and produces a dictionary containing
pairwise divergence and location of all SNPs for all possible sequence pairs.
'''

import argparse
from Bio import SeqIO
import numpy as np
from itertools import combinations
import pickle

def shared_strains(align_list):
    """
    Return list of strains shared across all genome segments.
    """
    setlist = []
    for align in align_list:
        strains = set()
        for record in SeqIO.parse(align, "fasta"):
            strains.add(record.id)
        setlist.append(strains)
    shared = set.intersection(*setlist)
    return shared

    return mapping, genome_length

def concat(shared, align_list):
    """
    Concatenates genome segments together
    """
    mapping = {}
    for strain in shared:
        mapping[strain] = ""
    for align in align_list:
        for record in SeqIO.parse(align, "fasta"):
            if record.id in mapping:
                seq = str(record.seq)
                mapping[record.id] += seq
    return mapping

def convert(genome_dict):
    """
    It replaces nucleotide bases with numbers in order to find SNPs efficiently computationally and deal with ambiguous nucleotide bases.
    C,T,G,A are replaced with the first 4 prime numbers. Ambiguous nucleotide bases are replaced with the multiple of the prime numbers
    for whichever bases it can represent. For example, R, which could be A or G, is 35 because A is 7, and G is 5.
    """
    mapping = {}
    for record in genome_dict:
        seq = list(genome_dict[record])
        seq_conversion = {
            "C": 2,
            "T": 3,
            "G": 5,
            "A": 7,
            "D": 105,
            "R": 35,
            "Y": 6,
            "H": 42,
            "K": 15,
            "V": 70,
            "N": 210,
            "M": 14,
            "W": 21,
            "S": 10
        }
        sequence = [seq_conversion.get(base) for base in seq]
        array = np.asarray(sequence, dtype=int)
        mapping[record] = array
    return mapping

def compare(array1, array2):
    condition1 = (array1 != array2)
    condition2 = (array1 % array2 != 0)
    condition3 = (array2 % array1 != 0)
    location = np.where(condition1 & condition2 & condition3)
    return location[0]

def pairwise(mapping):
    """
    Compares things pairwise
    """
    counter = 0
    interval = 5000
    pairwise_dict = {}
    strain_pairs = list(combinations(mapping.keys(), 2))
    for pair in strain_pairs:
        if counter % interval == 0:
            print("[", end="")
            for x in range(int(counter / interval)):
                print("-", end="")
            for x in range(int(len(strain_pairs) / interval) - int(counter / interval)):
                print(" ", end="")
            print("]")
        genome_length = len(mapping[pair[0]])
        location = compare(mapping[pair[0]], mapping[pair[1]])
        distance = np.size(location)
        pairwise_dict[pair] = {
            "snp_loc": location,
            "distance": distance,
            "length": genome_length
        }
        counter += 1
    return pairwise_dict

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Pairwise analysis for each sequence.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--alignment', type=str, nargs='+', required=True, help="list of alignments")
    parser.add_argument('--shared', type=str, required=True, help = "name of output of shared strains")
    parser.add_argument('--dictionary', type=str, required=True, help = "name of output dictionary")
    args = parser.parse_args()

    #Create list of strains shared across all segments
    shared = shared_strains(args.alignment)

    #Save shared strains as pickle files
    with open(args.shared, 'wb') as file:
        pickle.dump(shared, file)

    #Concatenate segments into one genome
    concatenated = concat(shared, args.alignment)

    #Converts bases from letters to numbers for comparison
    converted = convert(concatenated)

    #Calculates pairwise divergence for every strain
    comparison_dict = pairwise(converted)

    #Save dictionary as pickle file.
    with open(args.dictionary, 'wb') as file:
        pickle.dump(comparison_dict, file)
