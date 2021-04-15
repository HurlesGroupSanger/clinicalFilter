"""copyright"""

import logging

from filtering.inheritance_autosomal import AutosomalFilter
from filtering.inheritance_allosomal import AllosomalFilter
from filtering.compound_hets import CompoundHetScreen
from filtering.inheritance_report import InheritanceReport

class InheritanceFiltering(object):

    def __init__(self, variants_per_gene, family, genes, regions, trusted_variants):
        self.variants_per_gene = variants_per_gene
        self.family = family
        self.genes = genes
        self.regions = regions
        self.trusted_variants = trusted_variants
        self.parents = None

    def inheritance_filter(self):

        candidate_variants = {}
        candidate_variants['single_variants'] = {}
        candidate_variants['compound_hets'] = {}

        # inheritance report gives data to populate a matrix for every possible
        # combination of proband and parent GT and affected status for each mode of
        # inheritance. Some variants will be counted twice (genes which are both
        # mono and biallelic). Split into X and autosome
        # inheritance_report = create_blank_inheritance_report()
        inhreport = InheritanceReport()

        # lof_cqs = ['transcript_ablation', 'splice_donor_variant', 'stop_lost',
        #            'splice_acceptor_variant', 'stop_gained', 'frameshift_variant',
        #            'start_lost']

        if self.family.has_both_parents():
            self.parents = 'both'
        elif self.family.has_no_parents():
            self.parents = 'none'
        elif self.family.has_dad():
            self.parents = 'dad_only'
        elif self.family.has_mum():
            self.parents = 'mum_only'
        else:
            print(
                "Error: should not get here as family must either have both, one \
                 or no parents")
            logging.error("Can't calculate inheritance type - family error")
            exit(1)

        if self.genes:
            self.inheritance_filter_genes(candidate_variants, inhreport)

        if self.regions:
            # todo filtering for regions
            pass

        if self.trusted_variants:
            # todo filtering for trusted variants
            pass

        compoundhets = CompoundHetScreen(candidate_variants, self.family)
        screened_candidate_variants = compoundhets.screen_compound_hets()

        return screened_candidate_variants, inhreport.inheritance_report

    def inheritance_filter_genes(self, candidate_variants, Inhreport):
        # inheritance filters for use with a gene list
        for hgncid in self.variants_per_gene.keys():
            if hgncid in self.genes.keys():
                if self.genes[hgncid]['chr'] in ['X', 'Y']:
                    allosomalfiltering = AllosomalFilter(self, candidate_variants, Inhreport, hgncid)
                    allosomalfiltering.allosomal_filter()
                else:
                    autosomalfiltering = AutosomalFilter(self, candidate_variants, Inhreport, hgncid)
                    autosomalfiltering.autosomal_filter()
            else:
                for v in self.variants_per_gene[hgncid].keys():
                    logging.info(
                        v + " gene not in DDG2P: " +
                        self.variants_per_gene[hgncid][v]['child'].symbol)

