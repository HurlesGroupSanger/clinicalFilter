import logging

from utils.utils import common_elements
from utils.utils import add_compound_het_to_candidates
from utils.utils import add_single_var_to_candidates

class AllosomalFilter(object):

    def __init__(self, Inheritancefilter, candidate_variants, inheritance_report, hgncid):
        self.family = Inheritancefilter.family
        self.parents = Inheritancefilter.parents
        self.variants_per_gene = Inheritancefilter.variants_per_gene
        self.candidate_variants = candidate_variants
        self.inheritance_report = inheritance_report
        self.gene = Inheritancefilter.genes[hgncid]
        self.hgncid = hgncid

    def allosomal_filter(self):
        if self.parents == 'both':
            self.allosomal_both_parents()
        elif self.parents == 'none':
            self.allosomal_no_parents()
        else:
            self.allosomal_single_parent()

    def allosomal_both_parents(self):
        '''screens variants in X/Y where there are both parents and adds candidates
        to candidate_variants
        '''
        pass

    def allosomal_no_parents(self):
        '''screens variants in X/Y where there are no parents and adds candidates
        to candidate_variants
        '''
        pass

    def allosomal_single_parent(self):
        '''screens variants in X/Y where there is one parent and adds candidates
        to candidate_variants
        '''
        pass
