import gzip

def read_VCF(vcf):
    if vcf.endswith('gz'):
        with gzip.open(vcf, 'rt') as f:
            for line in f:
                yield line
    else:
        with open(vcf, 'r') as f:
            for line in f:
                yield line

def parse_header_VCF(vcf):
    VCF_header = {'contigs':{}, 'filters':[], 'header_len':0}
    for line in read_VCF(vcf):
        if line.startswith('##'):
            VCF_header['header_len'] += 1
            if line.startswith('##FILTER'):
                VCF_header['filters'].append(line.split('ID=')[1].split(',')[0])
            elif line.startswith('##contig'):
                contig = line.split('ID=')[1].split(',')[0]
                length = int(line.split('length=')[1].split('>')[0])
                VCF_header['contigs'][contig] = length
        elif line.startswith('#CHROM'):
            VCF_header['header_len'] += 1
            VCF_header['samples'] = line.strip().split('\t')[9:]
        else:
            break
    return VCF_header

