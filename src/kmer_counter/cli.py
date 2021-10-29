"""Module that contains the command line application."""

import gzip
import sys
import argparse
import py2bit
from kmer_counter.readers import *
from kmer_counter.counters import count_indels, count_non_indels

def main(args = None):
    """
    Run the main program.

    This function is executed when you type `kmer_counter` or `python -m kmer_counter`.

    Arguments:
        args: Arguments passed from the command line.

    Returns:
        An exit code.
    """
    parser = argparse.ArgumentParser(
        prog='kmer_counter',
        description='''
        Count k-mers at SNVs, indel breakpoints or in genomic regions
    ''')
    subparsers = parser.add_subparsers(dest='command', help='command. Must choose what kind of file you want to count.')
    parser.add_argument('-V', '--version', action='store_true')

    snv_parser = subparsers.add_parser('snv', 
        description='''Count k-mers at SNVs.
        The k-mer will be centered around the mutation base and if the reference base is G or T then the reverse_complement of the reference k-mer will be counted instead.
        '''
    )
    snv_parser.add_argument('ref_genome', help='Reference genome in 2bit format', type=str)
    snv_parser.add_argument('mutations', 
        type=argparse.FileType('r'), 
        help='A vcf-like file with SNVs. First four columns should be: Chrom, pos, ref, alt. '
        'Other columns are ignored. Non-SNV variants are ignored.')
    snv_parser.add_argument('-r', '--radius', type=int, metavar='R', default=1,
        help='Count the k-mer from R bases before a position to R bases '
             'after a position. For each position in the inputfile.')

    indel_parser = subparsers.add_parser('indel', description='Count k-mers at indels.')
    indel_parser.add_argument('ref_genome', help='Reference genome in 2bit format', type=str)
    indel_parser.add_argument('mutations',
        type=argparse.FileType('r'), 
        help='A sorted vcf-like file with indels. First four columns should be: Chrom, pos, ref, alt.'
        'Other columns are ignored. Non-indel variants are ignored. Indels should be left-aligned.')
    indel_parser.add_argument('type', choices=['ins', 'del_start', 'del_end', 'all', 'del'],
        help='What type of indel breakpoint do you want to count?')
    indel_parser.add_argument('-r', '--radius', default=1, type=int, 
        help='How many base pairs before indel_start_point or after indel_end_point should be included '
        'as context annotation.')
    indel_parser.add_argument('--sample', action="store_true", help='Randomly choose one of the possible positions instead of counting the expected (non integer) count for each possible position of an ambigously aligned indel.')
    #indel_parser.add_argument('-v', '--verbose', action='store_true')


    bg_parser = subparsers.add_parser('background', description='Count kmers in (regions of) a genome')
    bg_parser.add_argument('ref_genome',  type=str,
        help='Reference genome in 2bit format',)
    bg_parser.add_argument('--bed', type=str,
        help='bed-file describing regions that should be counted. May be gzipped.')
    #bg_parser.add_argument('--wig', type=str,
    #    help='wig-file describing regions that should be counted. May be gzipped. '
    #        'The context at a position will be weigthed by the value from the '
    #        'wig-file at that position. The output counts will thus be floats '
    #        'and not integers')
    bg_parser.add_argument('--all_autosomes', action="store_true",
        help='All parts of the autosomes will be counted')
    bg_parser.add_argument('-r', '--radius', type=int, metavar='R',
        help='Count the k-mer from R bases before a position to R bases '
             'after a position. For each position in the inputfile.')
    bg_parser.add_argument('--before_after', type=int, nargs=2, metavar=('X','Y'),
        help='count the k-mer from X bases before a position to Y bases after a position. '
        'For each position in the inputfile.')
    bg_parser.add_argument('--reverse_complement_method', type=str, choices=['none', 'middle', 'lexicographic', 'both'],
        help='"none" means that alle k-mers are counted unchanged. "middle" means that the reverse complement of a k-mer is counted if the middle position is not a "A" or "C". "lexicographic" means that the reverse_complement is counted if it has a smaller lexicographic order. Default is "middle" if --radius option is used and "lexicographic" if --before_after is used.')

    args = parser.parse_args(args)

    if args.version:
        from kmer_counter import __version__
        print("version:", __version__)
        print()
        return 0

    if args.command not in ['snv', 'indel', 'background']:
        print('Error: must specify command.')
        print()
        parser.print_help()
        return 0

    if 'ref_genome' not in args:
        print('Error: ref_genome (as 2bit file) must be specified')
        print()
        parser.print_help()
        return 0

    tb = py2bit.open(args.ref_genome)

    if args.command == 'indel':
        kmer_count = count_indels(args.mutations, tb, args.type, args.radius, args.sample)
    elif args.command == 'snv':
        dreader = PosReader(args.mutations, tb)
        kmer_count = count_non_indels(tb, dreader, args.radius, args.radius, 'middle')     
    elif args.command == 'background':
        if args.radius is None == args.before_after is None:
            raise Exception('Either the --radius or the --before_after option should be used (not both).')
        if not args.radius is None:
            assert args.radius>0
            before = args.radius
            after = args.radius
            if args.reverse_complement_method is None:
                args.reverse_complement_method = "middle"
        else:
            before, after = args.before_after
            assert before>=0
            assert after>=0
            if args.reverse_complement_method is None:
                args.reverse_complement_method = "none"        

        if not args.bed is None:
            if args.all_autosomes:
                raise Exception('Either --bed or--all_autosomes option should be used. Not both.')
            if args.bed.endswith('.bed.gz'):
                dreader = BedReader(gzip.open(args.bed, 'rt'))
            elif args.bed.endswith('.bed'):
                dreader = BedReader(open(args.bed))
            elif args.bed == '-':
                dreader = BedReader(sys.stdin)
            else:
                raise Exception('bed file should end with ".bed" or ".bed.gz"')
        elif args.all_autosomes:
            dreader = AllAutoReader(tb)
        else:
            raise Exception('Either --bed or--all_autosomes option should be used')
        kmer_count = count_non_indels(tb, dreader, before, after, args.reverse_complement_method)
    
    for x in kmer_count:
        print(x, kmer_count[x])
    
    return 0
