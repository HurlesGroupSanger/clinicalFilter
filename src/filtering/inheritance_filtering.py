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

import logging

from filtering.inheritance_autosomal import AutosomalFilter
from filtering.inheritance_allosomal import AllosomalFilter


class InheritanceFiltering(object):
    """
    Inheritance filtering of SNVs/Indels
    """

    def __init__(self, variants_per_gene, family, genes, regions, trusted_variants, candidate_variants, inhreport):
        self.variants_per_gene = variants_per_gene
        self.family = family
        self.genes = genes
        self.regions = regions
        self.trusted_variants = trusted_variants
        self.candidate_variants = candidate_variants
        self.inhreport = inhreport
        self.parents = None
        if self.family.has_both_parents():
            self.parents = "both"
        elif self.family.has_no_parents():
            self.parents = "none"
        elif self.family.has_dad():
            self.parents = "dad_only"
        elif self.family.has_mum():
            self.parents = "mum_only"

    def inheritance_filter(self):
        if self.genes:
            self.inheritance_filter_genes()
        if self.regions:
            # todo filtering for regions
            pass

        if self.trusted_variants:
            # todo filtering for trusted variants
            pass

    def inheritance_filter_genes(self):
        """
        Inheritance filters for use with a gene list
        """
        for hgncid in self.variants_per_gene.keys():
            if hgncid in self.genes.keys():
                if self.genes[hgncid]["chr"] in ["X", "Y"]:
                    allosomalfiltering = AllosomalFilter(self, hgncid)
                    allosomalfiltering.allosomal_filter()
                else:
                    autosomalfiltering = AutosomalFilter(self, hgncid)
                    autosomalfiltering.autosomal_filter()
            else:
                for v in self.variants_per_gene[hgncid].keys():
                    logging.info(v + " gene not in DDG2P: " + self.variants_per_gene[hgncid][v]["child"].symbol)
