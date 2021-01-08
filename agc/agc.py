#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Lapeyre Adrien"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Lapeyre Adrien"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Lapeyre Adrien"
__email__ = "lapeyreadr@eisti.eu"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default=400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default=10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default=100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default=8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file, minseqlen):
    """Returns a generator of sequences
      :Parameters:
          amplicon_file: Path to the amplicon_file
          minseqlen: Minimal length of sequences (int)
      Returns: generator of sequences
    """
    # Open the file in read mode
    with gzip.open(amplicon_file, 'rt') as my_file:
        lines = my_file.readlines()
        sequence = ""
        # Start on the second line
        for line in lines[1:]:
            # If this is the "title" of the sequence
            if line[0] == ">":
                # Yield the sequence if her length is enough
                if len(sequence) >= minseqlen:
                    yield sequence
                # Start a new sequence
                sequence = ""
            else:
                # Add the line to the current sequence without the \n character
                sequence += line[:-1]
        # Yield the last sequence
        if len(sequence) >= minseqlen:
            yield sequence


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """Returns a generator of unique sequences
      :Parameters:
          amplicon_file: Path to the amplicon_file
          minseqlen: Minimal length of sequences (int)
          mincount: Minimum counting (int)
      Returns: generator of unique sequences
    """
    dic = {}
    # Sequence generator
    read = read_fasta(amplicon_file, minseqlen)
    for seq in read:
        # Add each seq in the dictionnary if it doesn't exist or +1 to his value
        if seq not in dic:
            dic[seq] = 1
        else:
            dic[seq] += 1
    # Run the dict in descending order
    for key, value in sorted(dic.items(), key=lambda x: x[1], reverse=True):
        # Yield the sequence if her counting is enough
        if value >= mincount:
            yield [key, value]


def get_chunks(sequence, chunk_size):
    """Returns a list of non overlapping sub sequences of length chunk_size
      :Parameters:
          sequence: Sequence (string)
          chunk_size: Sub sequences length (int)
      Returns: list of sub sequences (list)
    """
    sub_seqs = []
    # Build 4 sub sequences
    for i in range(4):
        sub_seqs.append(sequence[i*chunk_size:(i+1)*chunk_size])
    return sub_seqs


def get_unique(ids):
    """Returns a dict_keys object (with unique elements)"""
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    """Returns a list with elements which are in both two lists"""
    return list(set(lst1) & set(lst2))


def cut_kmer(sequence, kmer_size):
    """Returns a generator of kmers
      :Parameters:
          sequence: Sequence (string)
          kmer_size: Size of kmers (int)
      Returns: generator of kmers
    """
    # Yield each kmer of a sequence
    for i in range(len(sequence)-kmer_size+1):
        yield sequence[i:i+kmer_size]


def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    """Returns a dictionnary of kmers with the list of sequence ids of the sequences where they are
      :Parameters:
          kmer_dict: Dictionnary of kmers (dict)
          sequence: Sequence (string)
          id_seq: Id of the sequence (int)
          kmer_size: Size of kmers (int)
      Returns: dictionnary of kmers (dict)
    """
    # Kmer generator
    read = cut_kmer(sequence, kmer_size)
    for kmer in read:
        # Add each kmer in the dictionnary if it doesn't exist or his id_seq to the list
        if kmer not in kmer_dict:
            kmer_dict[kmer] = [id_seq]
        else:
            kmer_dict[kmer].append(id_seq)
    return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size):
    """Returns a list with the sequence ids of the 8 most similar sequences with the candidate one
      :Parameters:
          kmer_dict: Dictionnary of kmers (dict)
          sequence: Sequence (string)
          kmer_size: Size of kmers (int)
      Returns: list of ids (list)
    """
    result = []
    id_seq_string = ""
    # Kmer generator
    read = cut_kmer(sequence, kmer_size)
    for kmer in read:
        for key, value in kmer_dict.items():
            # Find equal kmer to the one used
            if kmer == key:
                # Add all the sequence ids of this kmer to the big string of ids
                for id_seq in value:
                    id_seq_string += str(id_seq)
    # Count the 8 most similar sequences
    counter = Counter(id_seq_string).most_common(8)
    # Add the id of those sequences to the result
    for elem in counter:
        result.append(int(elem[0]))
    return result


def get_identity(alignment_list):
    """Returns the percentage of identity of two sequences
      :Parameters:
          alignment_list: List with the two sequences to compare (list)
      Returns: percentage of identity of two sequences (float)
    """
    length = len(alignment_list[0])
    nb_equal = 0
    for i in range(length):
        if alignment_list[0][i] == alignment_list[1][i]:
            nb_equal += 1
    return round((nb_equal/length)*100, 2)


def detect_chimera(perc_identity_matrix):
    """Returns True if the candidate sequence is a chimera (False if not)
      :Parameters:
          perc_identity_matrix: Matrix giving for each segment the percentage of identity betweeen
                                a candidate sequence and mate sequences (list)
      Returns: (bool)
    """
    result = False
    test = True
    percentages = []
    # Run the matrix and create a list with all percentages
    for lign in perc_identity_matrix:
        for elem in lign:
            percentages.append(elem)
    # Compute the standard deviation
    std = statistics.stdev(percentages)
    # Check if percentages of identity are not equal for each segment and each parent
    if (perc_identity_matrix[0][0] == perc_identity_matrix[1][0] == perc_identity_matrix[2][0] ==
            perc_identity_matrix[3][0] and perc_identity_matrix[0][1] == perc_identity_matrix[1][1]
            == perc_identity_matrix[2][1] == perc_identity_matrix[3][1]):
        test = False
    # If the check is ok and the std > 5 then it is a chimera
    if std > 5 and test:
        result = True
    return result


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """Returns a generator of non chimera sequences
      :Parameters:
          amplicon_file: Path to the amplicon_file
          minseqlen: Minimal length of sequences (int)
          mincount: Minimum counting (int)
          chunk_size: Sub sequences length (int)
          kmer_size: Size of kmers (int)
      Returns: generator of non chimera sequences
    """
    non_chimera_seq_list = []
    # Sequence generator
    read = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    cpt = 0
    for seq, value in read:
        # The two first sequences are considered as non chimera
        if cpt == 0 or cpt == 1:
            non_chimera_seq_list.append([seq, value])
        # Then we evaluate other ones
        else:
            id_seq = 0
            kmer_dict = {}
            # First we build the kmer_dict with all non chimera sequences
            for target_seq in non_chimera_seq_list:
                # Segments generator
                chunk_list = get_chunks(target_seq[0], chunk_size)
                for chunk in chunk_list:
                    kmer_dict = get_unique_kmer(kmer_dict, chunk, id_seq, kmer_size)
                id_seq += 1
            # Then we build a list of ids of mate non chimera sequences
            mate_seq_list_id = search_mates(kmer_dict, seq, kmer_size)
            perc_identity_matrix = 4*[[]]
            # Then we compute the matrix with the percentages of identity
            for mate in mate_seq_list_id:
                # Segments generators
                chunk_list1 = get_chunks(seq, chunk_size)
                chunk_list2 = get_chunks(non_chimera_seq_list[mate][0], chunk_size)
                for i in range(len(chunk_list1)):
                    # Make alignment between two segments
                    alignment_list = nw.global_align(chunk_list1[i], chunk_list2[i], gap_open=-1,
                                                     gap_extend=-1, matrix=os.path.abspath(
                                                         os.path.join(os.path.dirname(__file__),
                                                                      "MATCH")))
                    # Compute their identity
                    identity = get_identity(alignment_list)
                    perc_identity_matrix[i].append(identity)
            # Finally we decide if the candidate sequence is a chimera or not
            chimera = detect_chimera(perc_identity_matrix)
            # If it is not, we add it to the non chimera sequences list
            if not chimera:
                non_chimera_seq_list.append([seq, value])
        cpt += 1
    # Yield each non chimera sequence with her counting
    for elem in non_chimera_seq_list:
        yield elem


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def write_OTU(OTU_list, output_file):
    pass


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()


if __name__ == '__main__':
    main()
