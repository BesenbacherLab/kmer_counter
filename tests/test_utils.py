import py2bit
from kmer_counter import utils

def test_indels():
    tb = py2bit.open('tests/data/test.2bit')
    #tb = py2bit.open('/Users/sobe/Data/2bit/hg38.2bit')
    #seq ='CTCGCGCGTT'
    #for i in range(1000000, 400700100):
    #    if i % 1000000 == 0:
    #        print(i)
    #    if seq == tb.sequence("chr1", i, i+len(seq)):
    #        print(i, seq)

    #test insertion
    possible_pos = utils.get_possible_indel_pos("chr1:1012200-1012400", 106, "C", "CCT", tb ,1)
    assert(possible_pos == [(106, 'CT'), (107, 'TC'), (108, 'CT')])
    kmers = []
    for rpos, ralt in possible_pos: 
        c1, c2 = utils.get_indel_contexts("chr1:1012200-1012400", rpos, 2, tb)
        kmers.append(c1)
    assert(kmers == ['ACCT', 'CCTT', 'CTTA'])
    
    #test deletion start
    possible_pos = utils.get_possible_indel_pos("chr1:40884100-40884200", 43, "TCG", "T", tb, 1)
    assert(possible_pos == [(43, 'CG'), (44, 'GC'), (45, 'CG'), (46, 'GC'), (47, 'CG')])
    kmers=[]
    for rpos, ralt in possible_pos: 
        c1, c2 = utils.get_indel_contexts("chr1:40884100-40884200", rpos, 2, tb)
        kmers.append(c1)
    assert(kmers == ['CTCG', 'TCGC', 'CGCG', 'GCGC', 'CGCG'])
    
    #test deltion end 
    kmers=[]
    for rpos, ralt in possible_pos:
        c1, c2 = utils.get_indel_contexts("chr1:40884100-40884200",rpos+(len("TCG")-len("T")), 2, tb)
        kmers.append(c1)
    assert(kmers == ['CGCG', 'GCGC', 'CGCG', 'GCGT', 'CGTT'])
    kmers=[]

def test_reverse_complement():
    s = 'ACCCGTTGA'
    res = 'TCAACGGGT'
    assert(utils.reverse_complement(s)==res)