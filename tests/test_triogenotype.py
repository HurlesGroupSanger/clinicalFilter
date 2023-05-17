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

import unittest

from tests.test_utils import create_test_snv
from tests.test_utils import create_test_person
from tests.test_utils import create_test_family

from variants.trio_genotype import add_trio_genotypes


class TestTrioGenotypes(unittest.TestCase):
    def setUp(self):
        """create variants and people"""
        self.maxDiff = None
        homaltvardata = {
            "chrom": "1",
            "pos": "100000",
            "ref": "A",
            "alt": "G",
            "consequence": "cq",
            "ensg": "ensg",
            "symbol": "KMTD2",
            "feature": "feature",
            "canonical": "YES",
            "mane": "MANE",
            "hgnc_id": "123",
            "max_af": "0",
            "max_af_pops": ".",
            "ddd_af": "0",
            "revel": "1",
            "polyphen": ".",
            "hgvsc": ".",
            "hgvsp": ".",
            "sex": "XY",
            "dnm": ".",
            "gt": "1/1",
            "gq": "50",
            "ac_XX": "2",
            "ac_XY": "5",
            "nhomalt_XX": "379",
            "nhomalt_XY": "4",
            "ddd_father_af": ".",
        }
        hetvardata = {
            "chrom": "1",
            "pos": "100000",
            "ref": "A",
            "alt": "G",
            "consequence": "cq",
            "ensg": "ensg",
            "symbol": "KMTD2",
            "feature": "feature",
            "canonical": "YES",
            "mane": "MANE",
            "hgnc_id": "123",
            "max_af": "0",
            "max_af_pops": ".",
            "ddd_af": "0",
            "revel": "1",
            "polyphen": ".",
            "hgvsc": ".",
            "hgvsp": ".",
            "sex": "XY",
            "dnm": ".",
            "gt": "0/1",
            "gq": "50",
            "ac_XX": "2",
            "ac_XY": "5",
            "nhomalt_XX": "379",
            "nhomalt_XY": "4",
            "ddd_father_af": ".",
        }
        homrefvardata = {
            "chrom": "1",
            "pos": "100000",
            "ref": "A",
            "alt": "G",
            "consequence": "cq",
            "ensg": "ensg",
            "symbol": "KMTD2",
            "feature": "feature",
            "canonical": "YES",
            "mane": "MANE",
            "hgnc_id": "123",
            "max_af": "0",
            "max_af_pops": ".",
            "ddd_af": "0",
            "revel": "1",
            "polyphen": ".",
            "hgvsc": ".",
            "hgvsp": ".",
            "sex": "XY",
            "dnm": ".",
            "gt": "0/0",
            "gq": "50",
            "ac_XX": "2",
            "ac_XY": "5",
            "nhomalt_XX": "379",
            "nhomalt_XY": "4",
            "ddd_father_af": ".",
        }
        self.homaltvar = create_test_snv(homaltvardata)
        self.hetvar = create_test_snv(hetvardata)
        self.homrefvar = create_test_snv(homrefvardata)

        self.child = create_test_person("fam", "child_id", "dad_id", "mum_id", "M", "2", "/vcf/path")
        self.mum = create_test_person("fam", "mum_id", "0", "0", "F", "1", "/vcf/path")
        self.dad = create_test_person("fam", "dad_id", "0", "0", "M", "1", "/vcf/path")

    def test_triogenotype_both_parents(self):
        family = create_test_family(self.child, self.mum, self.dad)
        geno_to_var = {"2": self.homaltvar, "1": self.hetvar, "0": self.homrefvar}
        # test every iteration fo child, mum and dad genotypes
        for childv_g in ["1", "2"]:
            for mumv_g in geno_to_var.keys():
                for dadv_g in geno_to_var.keys():
                    childv = geno_to_var[childv_g]
                    mumv = geno_to_var[mumv_g]
                    dadv = geno_to_var[dadv_g]
                    triostring = childv_g + mumv_g + dadv_g
                    self.childvar = childv
                    variants = {
                        "child": {"1_100000_A_G": self.childvar},
                        "mum": {"1_100000_A_G": mumv},
                        "dad": {"1_100000_A_G": dadv},
                    }
                    add_trio_genotypes(family, variants)
                    self.assertEqual(self.childvar.triogenotype, triostring)

    def test_triogenotype_dad_only(self):
        # create family with no mother
        family = create_test_family(self.child, None, self.dad)
        geno_to_var = {"2": self.homaltvar, "1": self.hetvar, "0": self.homrefvar}

        for childv_g in ["1", "2"]:
            for dadv_g in geno_to_var.keys():
                childv = geno_to_var[childv_g]
                dadv = geno_to_var[dadv_g]
                triostring = childv_g + "NA" + dadv_g
                self.childvar = childv
                variants = {"child": {"1_100000_A_G": self.childvar}, "mum": {}, "dad": {"1_100000_A_G": dadv}}
                add_trio_genotypes(family, variants)
                self.assertEqual(self.childvar.triogenotype, triostring)

    def test_triogenotype_mum_only(self):
        # create family with no father
        family = create_test_family(self.child, self.mum, None)
        geno_to_var = {"2": self.homaltvar, "1": self.hetvar, "0": self.homrefvar}

        for childv_g in ["1", "2"]:
            for mumv_g in geno_to_var.keys():
                childv = geno_to_var[childv_g]
                mumv = geno_to_var[mumv_g]
                triostring = childv_g + mumv_g + "NA"
                self.childvar = childv
                variants = {"child": {"1_100000_A_G": self.childvar}, "mum": {"1_100000_A_G": mumv}, "dad": {}}
                add_trio_genotypes(family, variants)
                self.assertEqual(self.childvar.triogenotype, triostring)

    def test_triogenotype_no_parents(self):
        family = create_test_family(self.child, None, None)
        geno_to_var = {"2": self.homaltvar, "1": self.hetvar, "0": self.homrefvar}

        for childv_g in ["1", "2"]:
            childv = geno_to_var[childv_g]
            triostring = childv_g + "NANA"
            self.childvar = childv
            variants = {"child": {"1_100000_A_G": self.childvar}, "mum": {}, "dad": {}}
            add_trio_genotypes(family, variants)
            self.assertEqual(self.childvar.triogenotype, triostring)


if __name__ == "__main__":
    unittest.main()
