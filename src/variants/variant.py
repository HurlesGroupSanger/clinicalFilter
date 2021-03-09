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
        self.inheritance_type = None
        self.triogenotype = None

        self.set_genotype()
        self.standardise_chromosome()
        self.set_inheritance_type()
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

    def parse_hgnc_id(self):
        '''strip HGNC: from hgnc_id'''
        if self.hgnc_id.startswith('HGNC:'):
            self.hgnc_id = self.hgnc_id[5:]

    def set_inheritance_type(self):
        self.inheritance_type = 'autosomal'
        if self.chrom == 'X' and self.gender == 'M':
            self.inheritance_type = "XChrMale"
        elif self.chrom == 'X' and self.gender == 'F':
            self.inheritance_type = "XChrFemale"
        elif self.chrom == 'Y' and self.gender == 'M':
            self.inheritance_type = "YChrMale"
        elif self.chrom == 'Y' and self.gender == 'F':
            self.inheritance_type = "YChrFemale"




