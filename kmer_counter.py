import sys
import argparse
import random
import py2bit
import gzip
from collections import defaultdict

complement = str.maketrans('ATCGN', 'TAGCN')
def reverse_complement(s):
    return s.translate(complement)[::-1]

WrongNumberOfInputException = \
    Exception('Exactly one of the input options --bed, --pos or --all_autosomes should be used')

class PosReader:
    def __init__(self, f, tb):
        self.f = f
        self.tb = tb
    def __iter__(self):
        for line in self.f:
            L = line.split()
            chrom, pos, ref, alt = L[:4]
            if not (len(ref) == 1 and len(alt) == 1):
                print("Warning. Non-SNV variant ignored:", line.strip(), file=sys.stderr)
                continue
            if ref != self.tb.sequence(chrom, int(pos)-1, int(pos)):
                print("Warning. Reference allele dosn't match 2bit file:", line.strip(), file=sys.stderr)
                continue
            yield chrom, int(pos)-1, 1

class BedReader():
    def __init__(self,f):
        self.f = f
    def __iter__(self):
        for line in self.f:
            L = line.split()
            chrom, start, end = L[:3]
            for pos in range(int(start), int(end)):
                yield chrom, pos, 1

class AllAutoReader():
    def __init__(self, tb):
        self.tb = tb
        if any(x.startswith('chr') for x in tb.chroms()):
            prefix = 'chr'
        else:
            prefix = ''
        self.autosomes = [prefix + str(x) for x in range(1,23)]
    def __iter__(self):
        for chrom in self.autosomes:
            start = 0
            end = self.tb.chroms()[chrom]
            for pos in range(int(start), int(end)):
                yield chrom, pos, 1


def get_indel_contexts(chrom, pos, radius, tb):
    context = tb.sequence(chrom, pos-radius, pos+radius)
    #context = genome[chrom.encode('utf-8')][pos-args.radius:pos+args.radius].upper()
    if 'N' in context or len(context) < radius*2:
        raise Exception('bad_context')
    context2 = reverse_complement(context)
    return (context, context2)


def get_possible_indel_pos(chrom, pos, ref, alt, tb, buffer=50):
    #print(ref, alt)
    if buffer <= len(alt)+len(ref):
        return get_possible_indel_pos(chrom, pos, ref, alt, tb, buffer*2)
    # This method assumes that the indel is left aligned
    # long_seq og short_seq starts at the first position not shared between
    # alt og ref
    if len(ref) > len(alt): # This variant is a deletion
        assert(len(alt)==1)
        assert(ref[0] == alt)
        #print(tb.sequence(chrom, pos-1, pos+len(ref)-1))
        long_seq = tb.sequence(chrom, pos, pos+buffer) # Reference sequence
        #print(long_seq[:len(ref)-1], ref[1:])
        assert(long_seq[:len(ref)-1] == ref[1:])
        #long_seq = genome[pos:pos+buffer] # Reference sequence
        #short_seq = long_seq[:1] + long_seq[len(ref):] # Alternative sequence
        short_seq = long_seq[len(ref)-1:] # Alternative sequence
        i= 0
    elif len(ref) < len(alt): # This variant is an Insertion
        assert(len(ref)==1)
        assert(ref == alt[0])
        #alternative = genome[pos-1:pos] + alt + genome[pos+1:pos+buffer]
        short_seq = tb.sequence(chrom, pos, pos+buffer)
        #short_seq = genome[pos:pos+buffer] # Reference seqeunce
        long_seq =  alt[1:] + short_seq # Alternative seqeunce
        #long_seq = alt + genome[pos+1:pos+buffer] # Alternative seqeunce
        i = 0 
    else:
        assert(False)
    L = []
    len_diff = len(long_seq) - len(short_seq)
    # print(len_diff, len(alt))
    # print(i)
    # print("short long")
    # print(short_seq)
    # print(long_seq)
    # print("short start long start")
    # print(short_seq[:i])
    # print(long_seq[:i])
    # if not (short_seq[:i] == long_seq[:i]):
    #     print('PROBLEM start')
    #     return [(i+pos, long_seq[i:i+len_diff])]
    # print("short end long end")
    # print(short_seq[i:])
    # print(long_seq[i+len_diff:])
    # if not (short_seq[i:] == long_seq[i+len_diff:]):
    #     print('PROBLEM end')
    #     return [(i+pos, long_seq[i:i+len_diff])]
    assert(short_seq[:i] == long_seq[:i] and short_seq[i:] == long_seq[i+len_diff:])
    while True:     
        #print("prefixes = ", short_seq[:i], long_seq[:i])
        #print("suffix = ", short_seq[i:], long_seq[i+len_diff:])
        if short_seq[:i] == long_seq[:i] and short_seq[i:] == long_seq[i+len_diff:]:
            #L.append((i+pos, long_seq[i:i+len_diff], genome[i+pos-2:i+pos+2]))            
            L.append((i+pos, long_seq[i:i+len_diff]))            
            i+=1
            if i + len_diff >= buffer:
                return get_possible_indel_pos(chrom, pos, ref, alt, tb, buffer*2)
        else:
            return L

def test_indels():
    tb = py2bit.open('/home/besen/Data/2bit/hg38.2bit')
    #tb = py2bit.open('/Users/sobe/Data/2bit/hg38.2bit')
    #seq ='CTCGCGCGTT'
    #for i in range(1000000, 400700100):
    #    if i % 1000000 == 0:
    #        print(i)
    #    if seq == tb.sequence("chr1", i, i+len(seq)):
    #        print(i, seq)

    print(get_possible_indel_pos("chr1", 1012306, "C", "CCT", tb ,1))
    # [(1012306, 'CT'), (1012307, 'TC'), (1012308, 'CT')]
    kmers = []
    for rpos, ralt in get_possible_indel_pos("chr1", 1012306, "C", "CCT", tb ,1): 
        c1, c2 = get_indel_contexts("chr1", rpos, 2, tb)
        kmers.append(c1)
    # ['ACCT', 'CCTT', 'CTTA']
    print(kmers)
    kmers=[]
    # [(40884143, 'CG'), (40884144, 'GC'), (40884145, 'CG'), (40884146, 'GC'), (40884147, 'CG')]
    print(get_possible_indel_pos("chr1", 40884143, "TCG", "T", tb, 1))
    for rpos, ralt in get_possible_indel_pos("chr1", 40884143, "TCG", "T", tb, 1): 
        c1, c2 = get_indel_contexts("chr1", rpos, 2, tb)
        kmers.append(c1)
    print(kmers)
    #['CTCG', 'TCGC', 'CGCG', 'GCGC', 'CGCG']
    kmers=[]
    for rpos, ralt in get_possible_indel_pos("chr1", 40884143, "TCG", "T", tb, 1): 
        c1, c2 = get_indel_contexts("chr1",rpos+(len("TCG")-len("T")), 2, tb)
        kmers.append(c1)
    print(kmers)
    # ['CGCG', 'GCGC', 'CGCG', 'GCGT', 'CGTT']
    kmers=[]

def count_indels(args):
    if args.sample:
        kmer_count = defaultdict(int)
    else:        
        kmer_count = defaultdict(float)

    for line in args.mutations:
        chrom, pos, ref, alt = line.split()[:4]
        if args.verbose:
            print(chrom, pos, ref, alt)#, file=sys.stderr)
        pos = int(pos)
        if len(ref) > len(alt):
            if len(alt) != 1:
                if (ref[-(len(alt)-1):] != alt[1:]):
                    print("Warning. Complex variant ignored:", line.strip(), file=sys.stderr)
                    continue
                #print("Variant before", alt, ref)
                ref = ref[:-(len(alt)-1)]
                alt = alt[:1]
                #print("Variant after", alt, ref)
        elif len(alt)>len(ref):
            if len(ref) != 1:
                if(ref[1:] != alt[-(len(ref)-1):]):
                    print("Warning. Complex variant ignored:", line.strip(), file=sys.stderr)
                    continue
                alt = alt[:-(len(ref)-1)]
                ref = ref[:1]
        if alt[0] != ref[0]:
            print("Warning. Not left-aligned variant ignored:", line.strip())
            continue
        L = get_possible_indel_pos(chrom, pos, ref, alt, tb)
        if args.sample:
            L = [random.choice(L)]
        for rpos, ralt in L:
            if len(ref) < len(alt): # This variant is an Insertion
                try:
                    c1, c2 = get_indel_contexts(chrom, rpos, args.radius, tb)
                except:
                    continue
                if args.type in ['ins', 'all']:
                    if args.sample:
                        kmer_count[c1] += 1
                        kmer_count[c2] += 1
                    else:
                        kmer_count[c1] += 1.0/(len(L)*2)
                        kmer_count[c2] += 1.0/(len(L)*2)
            elif len(ref) > len(alt): # This variant is a Deletion
                try:
                    c1, c2 = get_indel_contexts(chrom, rpos, args.radius, tb)
                except:
                    continue
                if args.type in ['del_start', 'all', 'del']:
                    if args.sample:
                        kmer_count[c1] += 1
                    else:
                        kmer_count[c1] += 1.0/(len(L)*4)
                if args.type in ['del_end', 'all', 'del']:
                    if args.sample:
                        kmer_count[c2] += 1
                    else:
                        kmer_count[c2] += 1.0/(len(L)*4)
                try:
                    c1, c2 = get_indel_contexts(chrom, rpos+(len(ref)-len(alt)), args.radius, tb)
                except:
                    continue
                if args.type in ['del_start', 'all', 'del']:
                    if args.sample:
                        kmer_count[c2] += 1
                    else:
                        kmer_count[c2] += 1.0/(len(L)*4)
                if args.type in ['del_end', 'all', 'del']:
                    if args.sample:
                        kmer_count[c1] += 1
                    else:
                        kmer_count[c1] += 1.0/(len(L)*4)
    return kmer_count

def count_non_indels(tb, dreader, before, after, reverse_complement_method):
    kmer_count = defaultdict(int)
    for chrom, pos, scale in dreader:
        if pos<before or pos > (tb.chroms()[chrom]-after):
            continue
        try:
            context = tb.sequence(chrom, pos-before, pos+after+1)
            if 'N' in context:
                continue
        except KeyError:
            print('Unknown position:', chrom, pos, file=sys.stderr)
            continue
        if reverse_complement_method == "middle" and context[before] in ['T','G']:
            context = reverse_complement(context)
        elif reverse_complement_method == "lexicographic":
            context2 = reverse_complement(context)
            if context2 < context:
                context = context2
        kmer_count[context] += 1
    return kmer_count

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='''
    Count k-mers at SNVs, indel breakpoints or in a background file
    ''')
    subparsers = parser.add_subparsers(dest='command', help='command. Must choose what kind of file you want to count.')
    parser.add_argument('-t', '--test', action='store_true')
    parser.add_argument('-v', '--verbose', action='store_true')
    #parser.add_argument('-l', '--log_file', type=argparse.FileType('w'))
    
    snv_parser = subparsers.add_parser('snv', description='Count k-mers at SNVs.')
    snv_parser.add_argument('ref_genome', help='Reference genome in 2bit format', type=str)
    snv_parser.add_argument('mutations', 
        type=argparse.FileType('r'), 
        help='A vcf-like file with SNVs. First four columns should be: Chrom, pos, ref, alt. '
        'Other columns are ignored. Non-SNV variants are ignored.')
    snv_parser.add_argument('-r', '--radius', type=int, metavar='R',
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
    indel_parser.add_argument('-v', '--verbose', action='store_true')


    bg_parser = subparsers.add_parser('background', description='Count kmers in a genome')
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
    bg_parser.add_argument('--reverse_complement_method', type=str, choices=['none', 'middle', 'lexicographic'],
        help='"none" means that alle k-mers are counted unchanged. "middle" means that the reverse complement of a k-mer is counted if the middle position is not a "A" or "C". "lexicographic" means that the reverse_complement is counted if it has a smaller lexicographic order. Default is "middle" if --radius option is used and "lexicographic" if --before_after is used.')

    args = parser.parse_args()
    if args.test:
        test_indels()
        sys.exit()

    tb = py2bit.open(args.ref_genome)

    if args.command == 'indel':
        kmer_count = count_indels(args)
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
                raise WrongNumberOfInputException
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
            raise WrongNumberOfInputException
        kmer_count = count_non_indels(tb,dreader, before, after, args.reverse_complement_method)

    for x in kmer_count:
        print(x, kmer_count[x])
