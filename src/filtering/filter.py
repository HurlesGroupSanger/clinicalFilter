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

from file_loading.load_genes_and_regions import load_genes
from file_loading.load_vcfs import load_variants
from variants.trio_genotype import add_trio_genotypes
from filtering.preinheritance_filtering import PreInheritanceFiltering
from filtering.inheritance_filtering import InheritanceFiltering
from filtering.inheritance_cnv import CNVFiltering
from filtering.postinheritance_filter import PostInheritanceFiltering
from filtering.inheritance_report import InheritanceReport
from filtering.compound_hets import CompoundHetScreen


class Filter(object):
    """
    Class for filtering variants
    """

    def __init__(self, family, known_genes, known_regions,
                 trusted_variants, outdir):
        self.family = family
        self.known_genes = known_genes
        self.known_regions = known_regions
        self.trusted_variants = trusted_variants
        self.outdir = outdir
        self.candidate_variants = None
        self.candidate_variants = {'single_variants': {}, 'compound_hets': {}}
        self.inhreport = None
        self.inhreport = InheritanceReport()
        self.screened_candidate_variants = {}

    def filter_trio(self):
        """
        filter each trio
        """
        # if genes, regions or variants files are present we can create a list
        # of regions to load and therefore load fewer variants
        vcfregions = set()
        genes = None
        regions = None
        trusted_variants = None

        if self.known_genes:
            genes = load_genes(self.known_genes)

        if self.known_regions:
            # TODO add regions to the vcfregions set and populate regions variable
            pass

        if self.trusted_variants:
            # TODO add trusted variants locations to the vcfregions set and populate
            # trusted_regions variable
            pass

        if len(vcfregions) > 0:
            variants = load_variants(self.family, self.outdir, vcfregions)
        else:
            variants = load_variants(self.family, self.outdir)

        # add trio genotypes for each variant
        add_trio_genotypes(self.family, variants)

        # preinheritance filters
        preinheritancefilter = PreInheritanceFiltering(variants)
        variants_per_gene = preinheritancefilter.preinheritance_filter()

        # inheritance filters for SNVs
        inheritancefilter = InheritanceFiltering(variants_per_gene, self.family,
                                                 genes, regions,
                                                 trusted_variants,
                                                 self.candidate_variants,
                                                 self.inhreport)
        # candidate_variants, inheritance_report = inheritancefilter.inheritance_filter()
        inheritancefilter.inheritance_filter()

        # inheritance filters for CNVs
        cnvfilter = CNVFiltering(variants, self.family, genes, regions,
                                 trusted_variants, self.candidate_variants)
        cnvfilter.cnv_filter()

        # compound het screen
        compoundhets = CompoundHetScreen(self.candidate_variants, self.family)
        compoundhets.screen_compound_hets()

        # post inheritance filters
        postinheritancefilter = PostInheritanceFiltering(
            self.candidate_variants, self.family)
        filtered_candidate_variants = postinheritancefilter.postinheritance_filter()

        return filtered_candidate_variants, self.inhreport
