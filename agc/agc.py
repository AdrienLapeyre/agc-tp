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
          minseqlen: Minimal length of sequences
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
          minseqlen: Minimal length of sequences
          mincount: Minimum counting
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
    # Runs the dict in descending order
    for key, value in sorted(dic.items(), key=lambda x: x[1], reverse=True):
        # Yield the sequence if her counting is enough
        if value >= mincount:
            yield [key, value]


def get_chunks(sequence, chunk_size):
    pass


def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    return list(set(lst1) & set(lst2))


def cut_kmer(sequence, kmer_size):
    pass


def get_identity(alignment_list):
    pass


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass


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
