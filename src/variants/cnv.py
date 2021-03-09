"""
copyright
"""
from variants.variant import Variant

class CNV(Variant):
    """class for CNVs"""

    def set_genotype(self):
        pass

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