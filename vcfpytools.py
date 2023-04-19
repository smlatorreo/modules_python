import gzip
from operator import itemgetter


def parser(vcf_file):
    if vcf_file.endswith('gz'):
        with gzip.open(vcf_file, 'rb') as f:
            for line in f:
                yield line.decode()
    else:
        with open(vcf_file, 'r') as f:
            for line in f:
                yield line

def get_header(vcf_file):
    vcf = parser(vcf_file)
    for line in vcf:
        if line.startswith('#'):
            yield line
        else:
            break

def get_chr_lengths(vcf_file):
    header = get_header(vcf_file)
    lengths = {}
    for line in header:
        if line.startswith('##contig'):
            ID = line.split('ID=')[1].split(',')[0]
            l = line.split('length=')[1].split('>')[0]
            lengths[ID] = int(l)
    return(lengths)

def get_body(vcf_file):
    vcf = parser(vcf_file)
    for line in vcf:
        if line.startswith('#'):
            continue
        else:
            yield line

def get_samples(vcf_file):
    header = get_header(vcf_file)
    for line in header:
        if line.startswith('#CHROM'):
            return line.rstrip().split('\t')[9:]

def get_genotypes_hap(vcf_file, samples, binary = False, positions = False, filtered = True):
    '''
    Just for haploids
    '''
    indices = [get_samples(vcf_file).index(i) for i in samples]
    indices = [i+9 for i in indices]
    for line in get_body(vcf_file):
        filt = line.split('\t')[6]
        if filtered == True and filt != 'PASS':
            continue
        chrom = line.split('\t')[0]
        loc = int(line.split('\t')[1])
        pos = (chrom, loc)
        alleles = {'0':line.split('\t')[3], '1':line.split('\t')[4], '.':'.'}
        genotps = itemgetter(*indices) (line.rstrip().split('\t'))
        if len(samples) > 1: # For >= 2 samples
            if binary == False:
                yield {'Position':pos,'Alleles':(alleles['0'],alleles['1']),'Genotypes':[alleles[genotype[0]] for genotype in genotps]}
            else:
                yield {'Position':pos, 'Genotypes':[genotype[0] for genotype in genotps]}
        else: # For 1 sample
            if binary == False:
                yield {'Position':pos,'Alleles':(alleles['0'],alleles['1']),'Genotypes':alleles[genotps[0]]}
            else:
                yield {'Position':pos, 'Genotypes':genotps[0]}

def get_genotypes_dip(vcf_file, samples, binary = False, phase = '/', filtered = True):
    '''
    By deafult assumes unphased genotypes. Otherwise change phase to '|'
    '''
    indices = [get_samples(vcf_file).index(i) for i in samples]
    indices = [i+9 for i in indices]
    for line in get_body(vcf_file):
        filt = line.split('\t')[6]
        if filtered == True and filt != 'PASS':
            continue
        alleles = {'0':line.split('\t')[3], '1':line.split('\t')[4], '.':'.'}
        if len(alleles['1']) > 1: # Just biallelic SNPs
            continue
        genotps = [i.split(':')[0].split(phase) for i in itemgetter(*indices) (line.rstrip().split('\t'))]
        if len(samples) > 1: # More than 1 sample
            if binary == False:
                yield [phase.join((alleles[genotps[g][0]], alleles[genotps[g][1]])) for g in range(len(genotps))]
            else:
                yield [phase.join(g) for g in genotps]
        else: # For 1 sample
            if binary == False:
                phase.join((alleles[genotps[0][0]], alleles[genotps[0][1]]))
            else:
                yield phase.join(genotps[0])


