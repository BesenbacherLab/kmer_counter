"""
Entry-point module, in case you use `python -m kmer_counter`.
"""
import sys

from kmer_counter.cli import main

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
