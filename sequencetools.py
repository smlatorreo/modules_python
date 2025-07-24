import numpy as np
from random import shuffle
from math import log, sqrt

def readfasta(file):
    """Reads a .fasta file and returns a dictionary
    in which the keys are the sequences IDs and
    the values are the sequences
    :param file: fasta file
    :return: Dictionary
    """
    sequences = {}
    with open(file, 'r') as f:
        for line in f:
            sequence = ''
            if line.startswith('>'):
                nameseq = line.split('>')[1].strip()
            else:
                sequence = sequence + line.strip()
            sequences[nameseq] = sequence
    return sequences


def simulate_seq(length, GC_content=0.5, as_str=True):
    """Given a certain length and optionally a GC content
    value, this function returns a simulated DNA sequence
    :param length: Length of the sequence (int)
    :param GC_content: Value from 0.0 to 1.0 (float) (default = 0.5)
    :param as_str: Return as a string (True; default) of as a numpy array (False)
    :return: String or numpy array
    """
    pCG = GC_content / 2
    pAT = (1 - GC_content) / 2
    bases = ('A', 'T', 'C', 'G')
    probs = [pAT, pAT, pCG, pCG]   
    sequence = np.random.choice(bases, size=length, p=probs)
    if as_str == True:
        return sequence.view(f'U{length}').item()
    else:
        return sequence

def rev_comp(seq):
    """It returns the reverse-complement sequence of
    an input sequence
    :param seq: Input sequence (str)
    """
    complement = {'A':'T','T':'A','C':'G','G':'C','N':'N','.':'.','?':'?'}
    revcomp = ''.join([complement[n] for n in seq][::-1])
    return revcomp


def kimura2parameter(Seq1, Seq2):
    """Given two nucleotide sequences of the same length,
    this function returns their Kimura two-parameter distance
    :param Seq1: Sequence 1 (str) or (list)
    :param Seq2: Sequence 2 (str) or (list)
    :return: float
    """
    if len(Seq1) != len(Seq2):
        raise ValueError('Sequences need to be with the same length')
    bases = ('A','C','G','T')
    transitions = {('C','T'),('T','C'),('A','G'),('G','A')}
    Seq1 = [i.upper() for i in Seq1]
    Seq2 = [i.upper() for i in Seq2]
    N = 0
    counts = {'t':0,'v':0}
    for n1, n2 in zip(Seq1, Seq2):
        if n1 in bases and n2 in bases:
            N += 1
            if n1 != n2:
                if (n1, n2) in transitions:
                    counts['t'] += 1
                else:
                    counts['v'] += 1
    P = counts['t']/N
    Q = counts['v']/N
    try:
        d = (-1/2) * log((1 - (2*P) - Q) * sqrt(1 - (2*Q)))
    except ValueError:
        raise ValueError('K2P is not suitable for highly diverged sequences')
    return d
