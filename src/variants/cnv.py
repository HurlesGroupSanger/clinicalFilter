"""
Copyright (c) 2021 Genome Research Limited
Author: Ruth Eberhardt <re3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

from variants.variant import Variant


class CNV(Variant):
    """
    CNVs
    """

    def __init__(self, vardata):
        super().__init__(vardata)
        self.reportable_symbol = []
        self.reportable_hgnc_id = []

    def __repr__(self):
        # return str(self.__dict__)
        return (
            'CNV(chrom="{}", pos="{}", cnv_end="{}", ref="{}", alt="{}", '
            'cnv_type="{}", cnv_length="{}", cnv_filter="{}", '
            'consequence="{}", hgnc_id_all="{}", symbol_all="{}", '
            'gt="{}", cn="{}", cnv_inh="{}", sex="{}", genotype="{}", '
            'triogenotype="{}")'.format(
                self.chrom,
                self.pos,
                self.cnv_end,
                self.ref,
                self.alt,
                self.cnv_type,
                self.cnv_length,
                self.cnv_filter,
                self.consequence,
                self.hgnc_id_all,
                self.symbol_all,
                self.gt,
                self.cn,
                self.cnv_inh,
                self.sex,
                self.get_genotype(),
                self.triogenotype,
            )
        )

    def set_genotype(self):
        """
        Converts genotype to 0/1/2
        """
        pass

    def get_genotype(self):
        return self.genotype

    def set_triogenotype(self, triogeno):
        self.triogenotype = triogeno

    def is_snv(self):
        """
        Checks whether the variant is an SNV
        """
        return False

    def is_cnv(self):
        """
        Checks whether the variant is a CNV
        """
        return True

    def get_mum_genotype(self):
        if self.triogenotype[3:6] == "REF":
            mum_genotype = "0"
        else:
            mum_genotype = "1"
        return mum_genotype

    def get_dad_genotype(self):
        if self.triogenotype[6:9] == "REF":
            dad_genotype = "0"
        else:
            dad_genotype = "1"
        return dad_genotype

    def is_het(self):
        pass

    def is_hom_alt(self):
        pass

    def is_hom_ref(self):
        pass
