"""Tests for the `cli` module."""

import pytest

from kmer_counter import cli


def test_main():
    """Basic CLI test."""
    assert cli.main([]) == 0


def test_show_help(capsys):
    """
    Show help.

    Arguments:
        capsys: Pytest fixture to capture output.
    """
    with pytest.raises(SystemExit):
        cli.main(["-h"])
    captured = capsys.readouterr()
    assert "kmer_counter" in captured.out


def test_snv(capsys):
    res = {
        'GAT': 1,
        'CCC': 1,
        'CAA': 1,
        'TCC': 1,
        'AAA': 1,
        'GCG': 2,
        'CCA': 2,
        'GCA': 1,
        'ACC': 2,
        'CCG': 1
    }
    assert(cli.main(["snv","tests/data/test.2bit","tests/data/pos_test.txt","-r1"])==0)
    captured = capsys.readouterr()
    D = {}
    for line in captured.out.split('\n'):
        if line == "":
            continue
        kmer, count = line.split(" ")
        assert res[kmer] == int(count)    
    

def test_background(capsys):
    res = {
        'GCA': 7,
        'CAC': 16,
        'ACC': 15,
        'CCA': 14,
        'ACG': 3,
        'GCG': 12,
        'GCC': 24,
        'CCC': 14,
        'CAG': 11,
        'GCT': 7,
        'TAG': 8,
        'TAA': 6,
        'AAT': 5,
        'AAA': 8,
        'ACT': 4,
        'TAC': 1,
        'TCT': 4,
        'GAG': 9,
        'GAT': 12,
        'CAT': 1,
        'AAC': 5,
        'GAA': 1,
        'GAC': 3,
        'TCA': 12,
        'CAA': 8,
        'TCC': 9,
        'CCT': 16,
        'AAG': 5,
        'CCG': 16,
        'ACA': 3,
        'TAT': 1,    
    }
    assert(cli.main(["background","tests/data/test.2bit","--bed", "tests/data/bg_test.bed","-r1"])==0)
    captured = capsys.readouterr()
    D = {}
    for line in captured.out.split('\n'):
        if line == "":
            continue
        kmer, count = line.split(" ")
        assert res[kmer] == int(count)    


