import random
import sys
from collections import defaultdict
from kmer_counter.utils import *

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
            print("Warning. Not left-aligned variant ignored:", line.strip(), file=sys.stderr)
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
            print('Warning. Unknown position ignored:', chrom, pos, file=sys.stderr)
            continue
        if reverse_complement_method == "middle" and context[before] in ['T','G']:
            context = reverse_complement(context)
        elif reverse_complement_method == "lexicographic":
            context2 = reverse_complement(context)
            if context2 < context:
                context = context2
        elif reverse_complement_method == "both":
            context2 = reverse_complement(context)
            kmer_count[context2] += 1
        kmer_count[context] += 1
    return kmer_count