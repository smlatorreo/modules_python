from random import shuffle

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
    return(sequences)

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
    return(''.join(seq))

