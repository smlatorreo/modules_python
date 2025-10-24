class vcf_object:

    def __init__(self, path):
        self.path = path
        self.is_gzipped = True if path.endswith('gz') else False

    def dump(self):
        if self.is_gzipped == True:
            import gzip
            with gzip.open(self.path, 'rt') as f:
                for line in f:
                    yield line.strip()
        else:
            with open(self.path, 'r') as f:
                for line in f:
                    yield line

    def get_header(self):
        for line in self.dump():
            if line.startswith('#'):
                yield line
            else:
                break

    def get_contig_lengths(self):
        lengths = {}
        for line in self.get_header():
            if line.startswith('##contig'):
                ID = line.split('ID=')[1].split(',')[0]
                l = line.split('length=')[1].split('>')[0]
                lengths[ID] = int(l)
        return lengths

    def get_samples(self):
        for line in self.get_header():
            if line.startswith('#CHROM'):
                return line.split('\t')[9:]

    def get_body(self):
        for line in self.dump():
            if not line.startswith('#'):
                yield line

    def get_genotypes_hap(self, samples = [], binary = False, filtered = True):
        if len(samples) == 0:
            samples = self.get_samples()
        indices = [self.get_samples().index(i) + 9 for i in samples]
        for line in self.get_body():
            if filtered == True and line.split('\t')[6] != 'PASS':
                continue
            chrom = line.split('\t')[0]
            loc = int(line.split('\t')[1])
            pos = (chrom, loc)
            alleles = {'0':line.split('\t')[3], '1':line.split('\t')[4], '.':'.'}
            genotps = [line.split('\t')[i].split(':')[0] for i in indices]
            if binary == False:
                yield {'Position':pos, 'Alleles':(alleles['0'],alleles['1']), 'Genotypes':[alleles[g] for g in genotps]}
            else:
                yield {'Position':pos, 'Alleles':(alleles['0'],alleles['1']), 'Genotypes':[g for g in genotps]}

    def get_genotypes_dip(self, samples = [], phased = False, binary = False, filtered = True):
        if len(samples) == 0:
            samples = self.get_samples()
        phase = '/' if phased == False else '|'
        indices = [self.get_samples().index(i) + 9 for i in samples]
        for line in self.get_body():
            if filtered == True and line.split('\t')[6] != 'PASS':
                continue
            chrom = line.split('\t')[0]
            loc = int(line.split('\t')[1])
            pos = (chrom, loc)
            alleles = {'0':line.split('\t')[3], '1':line.split('\t')[4], '.':'.'}
            if len(alleles['1']) > 1: # Just biallelic SNPs
                continue
            genotps = [line.split('\t')[i].split(':')[0] for i in indices]
            if binary == False:
                nonbingeno = [phase.join([alleles[g.split(phase)[0]], alleles[g.split(phase)[1]]]) for g in genotps]
                yield {'Position':pos, 'Alleles':(alleles['0'],alleles['1']), 'Genotypes':nonbingeno}
            else:
                yield {'Position':pos, 'Alleles':(alleles['0'],alleles['1']), 'Genotypes':genotps}

    def get_info_values(self):
        info = [i for i in (self.get_header()) if i.startswith('##INFO')]
        metrics = [i.split('ID=')[1].split(',')[0] for i in info]
        def extract_metrics_info(info_col, metrics):
            metrics_out = []
            for metric in metrics:
                m = 'NA'
                if metric in info_col:
                    m = info_col.split(metric+'=')[1].split(';')[0]
                metrics_out.append(m)
            return metrics_out
        yield metrics
        for record in self.get_body():
            info_col = record.split('\t')[7]
            yield extract_metrics_info(info_col, metrics)

    def __str__(self):
        return 'VCF object\nFile name: \'{}\'\nGzipped: {}'.format(self.path, self.is_gzipped)


