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
pipx install kmer_counter
```

## Usage 

### Counting k-mers at SNVs
To count the 3-mers at SNVs do:
```bash
kmer_counter snv {genome}.2bit {snv_file}
```
Where the `{snv_file}` should be a vcf-like text file where the first four columns are: Chromosome, Position, Ref_Allele, Alt_Allele. Fx:

```
chr1  1000000  A G
chr1  1000200  G C
chr1  1000300  A T
chr1  1000500  C G
```
Comments or header lines starting with "#" are allowed and will be ignored and any additional columns are also allowed but ignored. So a vcf file is also a valid input file.
The Ref_Allele column should match the reference genome provided by the 2bit file. 2bit files can be downloaded from:
`https://hgdownload.cse.ucsc.edu/goldenpath/{genome}/bigZips/{genome}.2bit` where `{genome}` is a valid UCSC genome assembly name (fx. "hg38").

### Counting k-mers in genomic regions
To count all 5-mers in a bed file called {regions}.bed do:
```
kmer_counter background --bed {regions}.bed -radius 2 {genome}.2bit
```
By default all k-mers where the middle base is not A or C will be reverse complemented before being counted. This behaviour can be changed using the "--reverse_complement_method".
If we instead wants to count 4-mers, we can use the "--before_after" option: 
```
kmer_counter background --bed {regions}.bed --before_after 2 1 {genome}.2bit
```
When this option is used the default is not to reverse complement any of the k-mers but count all.

### Counting k-mers at indels
To count one of the possible insertion breakpoint 4-mer for each insertion in a vcf-like file with variants do:
```
kmer_counter indel -r 2 --sample {genome}.2bit {variants} ins
```
And for deletion start breakpoints:
```
kmer_counter indel -r 2 --sample {genome}.2bit {variants} del_start
```
This will produce 2 counts for each deletion; one for the start breakpoint and one for the reverse complement at the k-mer at end breakpoint.