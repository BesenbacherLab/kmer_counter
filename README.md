# kmer_counter

Count kmers in regions or at SNVs or at indel breakpoints.

## Requirements

kmer_counter requires Python 3.7 or above.

## Installation

With `pip`:
```bash
pip install kmer_counter
```

With [`pipx`](https://github.com/pipxproject/pipx):
```bash
pipx install --python python3.7 kmer_counter
```

## Usage 

### Counting k-mers at SNVs
To count the 3-mers at SNVs do:
```
kmer_counter snv {genome}.2bit {snv_file}
```
Where the `{snv_file}` should be a vcf-like text file where the first four columns are: Chromosome, Position, Ref_Allele, Alt_Allele. Fx:

```
chr1  1000000  A G
chr1  1000200  G C
chr1  1000300  A T
chr1  1000500  C G
```
Comments or headers lines starting with "#" are allowed and will be ignored and any additional columns are also allowed but ignored. So a vcf file is also a valid input file.
The Ref_Allele column should match the reference genome provided by the 2bit file. 2bit files can be downloaded from:
`https://hgdownload.cse.ucsc.edu/goldenpath/{genome}/bigZips/{genome}.2bit` where `{genome}` is a valid UCSC genome assembly name (fx. "hg38").

### Counting k-mers in genomic regions


### Counting k-mers at indels



