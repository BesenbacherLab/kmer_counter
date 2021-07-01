#!/bin/bash

#tmpfile1=$(gmktemp)
#echo chr1 1012200 1012400 | ~/Scripts/bed_add_seq.py ~/Data/2bit/hg38.2bit --fasta > $tmpfile1
#echo chr1 40884100 40884200 | ~/Scripts/bed_add_seq.py ~/Data/2bit/hg38.2bit --fasta >> $tmpfile1
#faToTwoBit $tmpfile1 tests/data/test.2bit

conda activate old
gseq 1012200 1012400 | sort -R |\
  head -20 | awk '{print "chr1",$1}' |\
  ~/Scripts/pos_add_sequence.py -r0 --two_bit ~/Data/2bit/hg38.2bit --reverse_complement_method none -u |\
  awk '$3!="T" {print "chr1:1012200-1012400", $2-1012200, $3, "T"}' | ~/Scripts/gorsort.sh > tests/data/pos_test.txt
cat tests/data/pos_test.txt | awk '{print $1,$2-10,$2+10}' > tests/data/bg_test.bed

#poetry run kmer_counter background tests/data/test.2bit --bed tests/data/bg_test.bed -r1 | gawk '{print "'\''"$1"'\'':",$2","}'
