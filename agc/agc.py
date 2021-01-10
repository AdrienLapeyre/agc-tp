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
    id_seq_string = ""
    # Kmer generator
    read = cut_kmer(sequence, kmer_size)
    for kmer in read:
        # Find equal kmer to the one used
        if kmer in kmer_dict:
            # Add all the sequence ids of this kmer to the big string of ids
            for id_seq in kmer_dict[kmer]:
                id_seq_string += str(id_seq)
    # Count the 8 most similar sequences
    counter = Counter(id_seq_string).most_common(8)
    # Return a list with their id
    return [int(elem[0]) for elem in counter]


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
    std = 0
    # Run the matrix and compute the mean standard deviation
    for lign in perc_identity_matrix:
        std += statistics.stdev(lign)
    std /= len(perc_identity_matrix)
    # Check if percentages of identity are not equal for each segment and each parent
    if not(perc_identity_matrix[0][0] == perc_identity_matrix[1][0] == perc_identity_matrix[2][0] ==
           perc_identity_matrix[3][0] and perc_identity_matrix[0][1] == perc_identity_matrix[1][1]
           == perc_identity_matrix[2][1] == perc_identity_matrix[3][1]):
        # If the std > 5 then it is a chimera
        if std > 5:
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
    id_seq = 0
    kmer_dict = {}
    # Sequence generator
    read = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    # Evaluate each sequence
    for seq, value in read:
        mate_seq_list_id = []
        # List of segments
        chunk_list = list(get_chunks(seq, chunk_size))
        # Build a list of ids of mate non chimera sequences for each segment
        for chunk in chunk_list:
            mate_seq_list_id.append(search_mates(kmer_dict, chunk, kmer_size))
        # Find parent sequences if there are
        parent_seq_list_id = common(common(mate_seq_list_id[0], mate_seq_list_id[1]),
                                    common(mate_seq_list_id[2], mate_seq_list_id[3]))
        perc_identity_matrix = [[], [], [], []]
        chimera = False
        # If there are at least 2 parents
        if len(parent_seq_list_id) >= 2:
            # Then we compute the matrix with the percentages of identity
            for parent in parent_seq_list_id[:2]:
                # List of segments of the parent
                chunk_list_p = list(get_chunks(non_chimera_seq_list[parent][0], chunk_size))
                for i in range(len(chunk_list)):
                    # Make alignment between two segments
                    alignment_list = nw.global_align(chunk_list[i], chunk_list_p[i], gap_open=-1,
                                                     gap_extend=-1, matrix=os.path.abspath(
                                                         os.path.join(os.path.dirname(__file__),
                                                                      "MATCH")))
                    # Compute their identity
                    identity = get_identity(alignment_list)
                    perc_identity_matrix[i].append(identity)
            # Finally we check if the candidate sequence is a chimera or not
            chimera = detect_chimera(perc_identity_matrix)
        # If it is not
        if not chimera:
            # We add it to the non chimera sequences list
            non_chimera_seq_list.append([seq, value])
            # We also add it to kmer_dict
            for chunk in chunk_list:
                kmer_dict = get_unique_kmer(kmer_dict, chunk, id_seq, kmer_size)
            id_seq += 1
            # And we yield it with her counting
            yield [seq, value]


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """Returns a list of OTU
      :Parameters:
          amplicon_file: Path to the amplicon_file
          minseqlen: Minimal length of sequences (int)
          mincount: Minimum counting (int)
          chunk_size: Sub sequences length (int)
          kmer_size: Size of kmers (int)
      Returns: list of OTU (list)
    """
    result = []
    # Non chimera sequences list (in ascending order)
    read = sorted(list(chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)),
                  key=lambda x: x[1])
    # Run the list (not the last sequence which is necessarily an OTU)
    cpt = 1
    for seq1, value1 in read[:-1]:
        otu = True
        # Run a second time to compare the sequence with all other (which have a bigger counting)
        for seq2, value2 in read[cpt:]:
            # Make alignment between two sequences
            alignment_list = nw.global_align(seq1, seq2, gap_open=-1, gap_extend=-1,
                                             matrix=os.path.abspath(os.path.join(
                                                 os.path.dirname(__file__), "MATCH")))
            # Check if the sequence is an OTU
            if seq1 != seq2 and get_identity(alignment_list) > 97 and value2 > value1:
                otu = False
                break
        # Add the sequence and her counting to the result if it is and OTU
        if otu:
            result.append([seq1, value1])
        cpt += 1
    # Add the last sequence to the list of OTU
    result.append([read[-1][0], read[-1][1]])
    # Sort the result in descending order
    result = sorted(result, key=lambda x: x[1], reverse=True)
    return result


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def write_OTU(otu_list, output_file):
    """Creates the OTU.fasta file
      :Parameters:
          otu_list: list of all OTU and their couting
          output_file: path of the OTU.fasta file
    """
    text = ''
    cpt = 1
    # Prepare the text to fasta format
    for elem in otu_list:
        text += ('>OTU_' + str(cpt) + ' occurrence:' + str(elem[1]) + '\n'
                 + fill(elem[0]) + '\n')
        cpt += 1
    # Write it in the output_file
    with open(output_file, 'w') as my_file:
        my_file.write(text)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # OTU search and creation of the OTU.fasta file
    otu_list = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount,
                                           args.chunk_size, args.kmer_size)
    write_OTU(otu_list, args.output_file)


if __name__ == '__main__':
    main()
