"""
copyright
"""
from variants.variant import Variant

class CNV(Variant):
    """class for CNVs"""

    def __init__(self, vardata):
        super().__init__(vardata)

    def __repr__(self):
        # return str(self.__dict__)
        return 'CNV(chrom="{}", pos="{}", cnv_end="{}", ref="{}", alt="{}", ' \
               'cnv_type="{}", cnv_length="{}", cnv_filter="{}", consequence="{}", ' \
               'hgnc_all="{}", symbol_all="{}", gt="{}", cn="{}", ' \
               'cnv_inh="{}", sex="{}", genotype="{}", triogenotype="{}")'.format(
            self.chrom, self.pos,
            self.cnv_end,
            self.ref, self.alt,
            self.cnv_type, self.cnv_length,
            self.cnv_filter, self.consequence,
            self.hgnc_id_all, self.symbol_all,
            self.gt, self.cn, self.cnv_inh,
            self.sex,
            self.get_genotype(),
            self.triogenotype)

    def set_genotype(self):
        '''converts genotype to 0/1/2'''
        pass
        # if self.cn == '0':
        #     self.genotype = '2'
        # elif self.cn == '1':
        #     self.genotype = '1'
        # elif self.cn == '3':
        #     self.genotype = '1'
        # elif int(self.cn) > 3:
        #     self.genotype = '2'
        # else:
        #     raise ValueError("invalid copy number")

    def get_genotype(self):
        return self.genotype

    def set_triogenotype(self, triogeno):
        self.triogenotype = triogeno

    def is_snv(self):
        """ checks whether the variant is for a CNV
        """
        return False

    def is_cnv(self):
        """ checks whether the variant is for a CNV
        """
        return True

    def is_het(self):
        pass

    def is_hom_alt(self):
        pass

    def is_hom_ref(self):
        pass