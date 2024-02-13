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
import logging


class SNV(Variant):
    """
    SNVs
    """

    def __init__(self, vardata):
        super().__init__(vardata)
        self.standardise_gt()
        self.calculate_ac_het_hemi()

    def __repr__(self):
        return (
            'SNV(chrom="{}", pos="{}", ref="{}", alt="{}", '
            'consequence="{}", protein_position="{}", ensg="{}", '
            'symbol="{}", feature="{}", canonical="{}", mane="{}", '
            'hgnc_id="{}", max_af="{}", max_af_pops="{}", ddd_af="{}", '
            'ddd_father_af="{}", revel="{}", polyphen="{}", hgvsc="{}", '
            'hgvsp="{}", dnm="{}", gt="{}", gq="{}", pid="{}", ad={}, '
            'sex="{}", genotype="{}", triogenotype="{}")'.format(
                self.chrom,
                self.pos,
                self.ref,
                self.alt,
                self.consequence,
                self.protein_position,
                self.ensg,
                self.symbol,
                self.feature,
                self.canonical,
                self.mane,
                self.hgnc_id,
                self.max_af,
                self.max_af_pops,
                self.ddd_af,
                self.ddd_father_af,
                self.revel,
                self.polyphen,
                self.hgvsc,
                self.hgvsp,
                self.dnm,
                self.gt,
                self.gq,
                self.pid,
                self.ad,
                self.sex,
                self.get_genotype(),
                self.triogenotype,
            )
        )

    def calculate_ac_het_hemi(self):
        """
        Calculate AC_het and AC_hemi
        """
        AC_hemi = 0
        AC_het = 0
        # change all '.' to '0'
        if self.ac_XX == ".":
            self.ac_XX = "0"
        if self.ac_XY == ".":
            self.ac_XY = "0"
        if self.nhomalt_XX == ".":
            self.nhomalt_XX = "0"
        if self.nhomalt_XY == ".":
            self.nhomalt_XY = "0"

        if self.chrom == "X" or self.chrom == "Y":
            AC_hemi = int(self.nhomalt_XY)
        else:
            total_AC = int(self.ac_XX) + int(self.ac_XY)
            total_nhom = int(self.nhomalt_XX) + int(self.nhomalt_XY)
            AC_het = total_AC - (2 * total_nhom)

        if AC_het < 0 or AC_hemi < 0:
            logging.info(
                "Warning - negative AC_het or AC_hemi for "
                + self.chrom
                + " "
                + self.pos
                + " AC_het="
                + str(AC_het)
                + " AC_hemi="
                + str(AC_hemi)
            )

        self.AC_het = str(AC_het)
        self.AC_hemi = str(AC_hemi)
        self.AC_gnomad = str(int(self.nhomalt_XX) + int(self.nhomalt_XY))

    def standardise_gt(self):
        """
        Reformat gt to ensure that lowest number allele is first and that
        the separator is /
        """
        newgt = self.gt
        newgt = newgt.replace("|", "/")
        gtsplit = list(newgt)
        if int(gtsplit[0]) > int(gtsplit[2]):
            newgt = gtsplit[2] + "/" + gtsplit[0]
        self.gt = newgt

    def set_genotype(self):
        """
        Converts genotype to 0/1/2
        """
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
        if genotype == "0":
            return True
        else:
            return False

    def is_dad_hom_ref(self):
        genotype = self.get_dad_genotype()
        if genotype == "0":
            return True
        else:
            return False

    def is_mum_hom_alt(self):
        genotype = self.get_mum_genotype()
        if genotype == "2":
            return True
        else:
            return False

    def is_dad_hom_alt(self):
        genotype = self.get_dad_genotype()
        if genotype == "2":
            return True
        else:
            return False

    def is_mum_het(self):
        genotype = self.get_mum_genotype()
        if genotype == "1":
            return True
        else:
            return False

    def is_dad_het(self):
        genotype = self.get_dad_genotype()
        if genotype == "1":
            return True
        else:
            return False

    def is_cnv(self):
        """
        Checks whether the variant is a CNV
        """
        return False

    def is_snv(self):
        """
        Checks whether the variant is an SNV
        """
        return True

    def is_het(self):
        """
        Is the variant a het?
        """
        if self.genotype == 1:
            return True
        else:
            return False

    def is_hom_alt(self):
        """
        Is the variant hom alt?
        """
        if self.genotype == 2:
            return True
        else:
            return False

    def is_hom_ref(self):
        """
        Is the variant hom ref?
        """
        if self.genotype == 0:
            return True
        else:
            return False
