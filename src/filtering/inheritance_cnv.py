import logging

class CNVFiltering(object):

    def __init__(self, variants, family, genes, regions, trusted_variants, candidate_variants):
        self.variants = variants
        self.family = family
        self.genes = genes
        self.regions = regions
        self.trusterd_variants = trusted_variants
        self.candidate_variants = candidate_variants
        self.parents = None
        if self.family.has_both_parents():
            self.parents = 'both'
        elif self.family.has_no_parents():
            self.parents = 'none'
        elif self.family.has_dad():
            self.parents = 'dad_only'
        elif self.family.has_mum():
            self.parents = 'mum_only'

    def cnv_filter(self):
        # filter CNVs here
        if self.genes:
            self.cnv_filter_genes()
        if self.regions:
            # todo filtering for regions
            pass
        if self.trusted_variants:
            # todo filtering for trusted variants
            pass

    def cnv_filter_genes(self):
        if self.parents == 'both':
            self.cnv_filter_parents()
        elif self.parents == 'none':
            self.cnv_filter_no_parents()
        else:
            self.cnv_filter_single_parent()

    def cnv_filter_parents(self):
        pass

    def cnv_filter_no_parents(self):
        pass

    def cnv_filter_single_parent(self):
        pass

    def cnv_inheritance_filter(self):
        pass

    def cnv_non_ddg2p_filter(self):
        pass

    def cnv_ddg2p_filter(self):
        pass