def parser(vcf_file):
    if vcf_file.endswith('gz'):
        import gzip
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

def get_genotypes(vcf_file, samples):
    from operator import itemgetter
    indices = [get_samples(vcf_file).index(i) for i in samples]
    indices = [i+9 for i in indices]
    for line in get_body(vcf_file):
        genotps = itemgetter(*indices) (line.rstrip().split('\t'))
        if len(samples) > 1: # For >= 2 samples
            yield [genotype[0] for genotype in genotps]
        else:
            yield genotps[0] # For 1 sample




