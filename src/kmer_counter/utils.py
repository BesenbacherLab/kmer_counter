complement = str.maketrans('ATCGN', 'TAGCN')
def reverse_complement(s):
    return s.translate(complement)[::-1]

def get_indel_contexts(chrom, pos, radius, tb):
    context = tb.sequence(chrom, pos-radius, pos+radius)
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
        long_seq = tb.sequence(chrom, pos, pos+buffer) # Reference sequence
        assert(long_seq[:len(ref)-1] == ref[1:])
        short_seq = long_seq[len(ref)-1:] # Alternative sequence
        i= 0
    elif len(ref) < len(alt): # This variant is an Insertion
        assert(len(ref)==1)
        assert(ref == alt[0])
        short_seq = tb.sequence(chrom, pos, pos+buffer) # Reference seqeunce
        long_seq =  alt[1:] + short_seq # Alternative seqeunce
        i = 0 
    else:
        assert(False)
    L = []
    len_diff = len(long_seq) - len(short_seq)
    assert(short_seq[:i] == long_seq[:i] and short_seq[i:] == long_seq[i+len_diff:])
    while True:     
        #print("prefixes = ", short_seq[:i], long_seq[:i])
        #print("suffix = ", short_seq[i:], long_seq[i+len_diff:])
        if short_seq[:i] == long_seq[:i] and short_seq[i:] == long_seq[i+len_diff:]:
            L.append((i+pos, long_seq[i:i+len_diff]))            
            i+=1
            if i + len_diff >= buffer:
                return get_possible_indel_pos(chrom, pos, ref, alt, tb, buffer*2)
        else:
            return L