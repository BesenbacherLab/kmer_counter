import sys
complement = str.maketrans('ATCGN', 'TAGCN')

def reverse_complement(s):
    """Reverse complement

    Args:
        s (str): input string should contain only A,C,G or T 

    Returns:
        str: reverse complement of input
    """
    return s.translate(complement)[::-1]

def get_indel_contexts(chrom, pos, radius, tb):
    """get sequence context around an indel breakpoint

    Args:
        chrom (str): chromosome
        pos (int): position
        radius (int): flank 
        tb: TwoBit object

    Raises:
        Exception: If context contains 'N' or length is smaller than expected

    Returns:
        (str,str): context and reverse complement of context.
    """
    context = tb.sequence(chrom, pos-radius, pos+radius)
    if 'N' in context or len(context) < radius*2:
        raise Exception('bad_context')
    context2 = reverse_complement(context)
    return (context, context2)


def get_possible_indel_pos(chrom, pos, ref, alt, tb, buffer=50):
    """Get list with possible positions of an indel
    This method assumes that the indel is left aligned.
    But does not check if that is true.

    Args:
        chrom (str): chromosome
        pos (int): position
        ref (str): reference base
        alt (str): alternative base
        tb (): TwoBit object
        buffer (int, optional): Number of upstream bases to read. Defaults to 50.
            If it turns out to be to small it will automatically be increased.

    Returns:
        list of int: all possible positions of the indel
    """

    if buffer <= len(alt)+len(ref):
        return get_possible_indel_pos(chrom, pos, ref, alt, tb, buffer*2)
    
    # long_seq og short_seq starts at the first position not shared between
    # alt og ref
    if len(ref) > len(alt): # This variant is a deletion
        assert(len(alt)==1)
        assert(ref[0] == alt)
        #print(chrom, pos, ref, alt)        
        long_seq = tb.sequence(chrom, pos, pos+buffer) # Reference sequence
        #print(long_seq)
        #print(long_seq[:len(ref)-1],ref[1:])
        #print(tb.sequence(chrom, pos-1, pos))
        if(long_seq[:len(ref)-1] != ref[1:]):
            print(f"WARNING. Reference does not match: {chrom} {pos} {ref} {alt}", file=sys.stderr)
            print(f"ref according to mutations file: {ref}", file=sys.stderr)
            oref = tb.sequence(chrom, pos-1, pos+len(ref)-1)
            print(f"ref according to 2bit file: {oref}", file=sys.stderr)
            print(f"Skipping variant.", file=sys.stderr)
            return None
        #assert(long_seq[:len(ref)-1] == ref[1:])
        short_seq = long_seq[len(ref)-1:] # Alternative sequence
        i = 0
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