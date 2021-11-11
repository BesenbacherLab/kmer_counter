"""
kmer_counter package.

Count kmers in regions or at SNVs or at indel breakpoints.
"""
from kmer_counter.utils import *

__all__: [reverse_complement, get_indel_contexts, get_possible_indel_pos]

__version__ = "0.2.1"
