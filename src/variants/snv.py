"""
copyright
"""
from variants.variant import Variant


class SNV(Variant):
    """class for SNVs"""

    def __init__(self, vardata):
        super().__init__(vardata)
        self.standardise_gt()

    def __repr__(self):
        # return str(self.__dict__)
        return 'SNV(chrom="{}", pos="{}", ref="{}", alt="{}", consequence="{}", ' \
               'protein_position="{}", ensg="{}", symbol="{}", feature="{}", canonical="{}", mane="{}", ' \
               'hgnc_id="{}", max_af="{}", max_af_pops="{}", ddd_af="{}", ' \
               'ddd_father_af="{}", revel="{}", polyphen="{}", hgvsc="{}", hgvsp="{}", dnm="{}", ' \
               'gt="{}", gq="{}", pid="{}", ad={}, sex="{}", ' \
               'genotype="{}", triogenotype="{}")'.format(
            self.chrom, self.pos,
            self.ref, self.alt,
            self.consequence,
            self.protein_position,
            self.ensg, self.symbol,
            self.feature,
            self.canonical, self.mane,
            self.hgnc_id, self.max_af,
            self.max_af_pops,
            self.ddd_af, self.ddd_father_af, self.revel,
            self.polyphen, self.hgvsc,
            self.hgvsp, self.dnm,
            self.gt, self.gq, self.pid, self.ad,
            self.sex,
            self.get_genotype(),
            self.triogenotype)

    def standardise_gt(self):
        '''Reformat gt to ensure that lowest number allele is first and that
        the separator is / '''
        newgt = self.gt
        newgt = newgt.replace('|', '/')
        gtsplit = list(newgt)
        if int(gtsplit[0]) > int(gtsplit[2]):
            newgt = gtsplit[2] + "/" + gtsplit[0]
        self.gt = newgt

    def set_genotype(self):
        '''converts genotype to 0/1/2'''
        if len(self.gt) != 3:
            raise ValueError("genotype should be three characters")
        else:
            gtsplit = list(self.gt)
            if gtsplit[0] == "0" and gtsplit[2] == "0":
                self.genotype = "0"
            elif gtsplit[0] == gtsplit[2]:
                self.genotype = "2"
            else:
                self.genotype = "1"

    def get_genotype(self):
        return self.genotype

    def set_triogenotype(self, triogeno):
        self.triogenotype = triogeno

    def get_mum_genotype(self):
        mum_genotype = self.triogenotype[1]
        return mum_genotype

    def get_dad_genotype(self):
        dad_genotype = self.triogenotype[2]
        return dad_genotype

    def is_mum_hom_ref(self):
        genotype = self.get_mum_genotype()
        if genotype == '0':
            return True
        else:
            return False

    def is_dad_hom_ref(self):
        genotype = self.get_dad_genotype()
        if genotype == '0':
            return True
        else:
            return False

    def is_mum_hom_alt(self):
        genotype = self.get_mum_genotype()
        if genotype == '2':
            return True
        else:
            return False

    def is_dad_hom_alt(self):
        genotype = self.get_dad_genotype()
        if genotype == '2':
            return True
        else:
            return False

    def is_mum_het(self):
        genotype = self.get_mum_genotype()
        if genotype == '1':
            return True
        else:
            return False

    def is_dad_het(self):
        genotype = self.get_dad_genotype()
        if genotype == '1':
            return True
        else:
            return False


    def is_cnv(self):
        """ checks whether the variant is for a CNV
        """
        return False

    def is_snv(self):
        """ checks whether the variant is for a CNV
        """
        return True


    def is_het(self):
        """is the variant a het?"""
        if self.genotype == 1:
            return True
        else:
            return False


    def is_hom_alt(self):
        """is the variant hom alt?"""
        if self.genotype == 2:
            return True
        else:
            return False


    def is_hom_ref(self):
        """is the variant hom ref?"""
        if self.genotype == 0:
            return True
        else:
            return False
