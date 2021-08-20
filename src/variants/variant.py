"""
copyright
"""
class Variant(object):
    """generic variant class, inherited by more specific classes such as CNV
     and SNV"""

    def __init__(self, vardata):
        for key in vardata:
            setattr(self, key, vardata[key])

        self.genotype = None
        self.triogenotype = None
        self.set_genotype()
        self.standardise_chromosome()
        self.parse_hgnc_id()

    def __eq__(self, other):
        return self.chrom == other.chrom and \
               self.pos == other.pos and \
               self.ref == other.ref and \
               self.alt == other.alt

    def standardise_chromosome(self):
        '''ensure chromosome is 1-22,X,Y'''
        if self.chrom.startswith('Chr') or self.chrom.startswith('chr'):
            self.chrom = self.chrom[3:]

    # def standardise_gt(self):
    #     '''standardise format of gt to always use / and have lowest number
    #     allele first'''


    def parse_hgnc_id(self):
        '''strip HGNC: from hgnc_id'''
        if self.hgnc_id.startswith('HGNC:'):
            self.hgnc_id = self.hgnc_id[5:]





