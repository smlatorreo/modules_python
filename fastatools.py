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

def simulate_seq(length, CG_content = 0.5):
    """Given a certain length and optionally a CG content
    value, this function returns a simulated DNA sequence
    :param length: Length of the sequence (int)
    :param CG_content: Value from 0.0 to 1.0 (float) (default = 0.5)
    :return: String
    """
    nucl = ('A','T','C','G')
    props = {'C':CG_content/2,'G':CG_content/2,'A':(1 -CG_content)/2,'T':(1 -CG_content)/2}
    seq = list(''.join([i*(int(length*props[i])) for i in nucl]))
    shuffle(seq)
    return ''.join(seq)

def kimura2parameter(Seq1, Seq2):
    if len(Seq1) != len(Seq2):
        print('Sequences need to be with the same length')
        return None
    bases = ('A','C','G','T')
    transitions = {('C','T'),('T','C'),('A','G'),('G','A')}
    Seq1 = [i.upper() for i in Seq1]
    Seq2 = [i.upper() for i in Seq2]
    counts = {'t':0,'v':0}
    for n1, n2 in zip(Seq1, Seq2):
        if n1 in bases and n2 in bases:
            if n1 != n2:
                if (n1, n2) in transitions:
                    counts['t'] += 1
                else:
                    counts['v'] += 1
    N = len(Seq1)
    P = counts['t']/N
    Q = counts['v']/N
    d = -0.5 * log((1 - 2 * P - Q) * sqrt(1 - 2 * Q))
    return d
