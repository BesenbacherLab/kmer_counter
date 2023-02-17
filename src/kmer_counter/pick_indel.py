import random
import sys
from collections import defaultdict
from kmer_counter.utils import *


def pick_random_indels(mutations, tb, radius):
    """ Pick a random indel position out of the possible indel positions.

    Args:
        mutations: file with mutations
        tb: TwoBit file.

    """    
    for line in mutations:
        chrom, pos, ref, alt = line.split()[:4]
        #if args.verbose:
        #    print(chrom, pos, ref, alt)#, file=sys.stderr)
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
        if L is None:
            continue
        
        rpos, ralt = random.choice(L)
        try:
            #new_ref =  tb.sequence(chrom, rpos, rpos+len(ref))
            #new_alt = tb.sequence(chrom, rpos, rpos+len(ref))
            c1, c2 = get_indel_contexts(chrom, rpos, radius, tb)
        except:
            continue
        new_ref = tb.sequence(chrom, rpos-1, rpos+len(ref)-1)
        if len(ref)<len(alt):
            new_alt = new_ref + ralt
            
        else:
            new_alt = new_ref[0]
        #assert(len(alt) == len(ralt)+1)
        #print(len(ref), len(alt), len(ralt))
        c = random.choice([c1,c2])
        print(chrom, rpos, new_ref, new_alt, c)

    return 0