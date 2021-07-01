import sys

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