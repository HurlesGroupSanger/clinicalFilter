import logging

from utils.utils import common_elements
from utils.utils import add_compound_het_to_candidates
from utils.utils import add_single_var_to_candidates


def allosomal_both_parents(hgncid, gene, variants, family, candidate_variants, inheritance_report,
                           lof_cqs):
    '''screens variants in X/Y where there are both parents and adds candidates
    to candidate_variants
    args:
    gene - info about the gene
    variants - all variants in this gene which pass preinheritance filters
    candidate_variants - output dict
    '''
    pass


def allosomal_no_parents(hgncid, gene, variants, candidate_variants, inheritance_report, lof_cqs):
    '''screens variants in X/Y where there are no parents and adds candidates
    to candidate_variants
    args:
    gene - info about the gene
    variants - all variants in this gene which pass preinheritance filters
    candidate_variants - output dict
    '''
    pass


def allosomal_single_parent(hgncid, gene, variants_per_gene, family,
                            candidate_variants, inheritance_report, lof_cqs):
    '''screens variants in X/Y where there is one parent and adds candidates
    to candidate_variants
    args:
    gene - info about the gene
    variants - all variants in this gene which pass preinheritance filters
    candidate_variants - output dict
    '''
    pass
