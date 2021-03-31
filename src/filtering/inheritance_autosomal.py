import logging

from utils.utils import common_elements
from utils.utils import add_compound_het_to_candidates
from utils.utils import add_single_var_to_candidates

def autosomal_both_parents(hgncid, gene, variants, family, candidate_variants, inheritance_report,
                           lof_cqs):
    '''screens variants in autosomes where there are both parents and adds
    candidates to candidate_variants
    args:
    gene - info about the gene
    variants - all variants in this gene which pass preinheritance filters
    candidate_variants - output dict
    '''
    for v in variants:
        # print(variants[v]['child'].genotype)
        # print(variants[v]['child'].gt)
        if variants[v]['child'].genotype == '1':
            #heterozygous
            for inh in gene['mode']:
                if inh == 'Biallelic':
                    biallelic_heterozygous_parents_filter(v, variants[v]['child'],
                                                     hgncid, inh, candidate_variants, inheritance_report)
                elif inh == 'Monoallelic':
                    monoallelic_heterozygous_parents_filter(v,
                                                          variants[v]['child'],
                                                          hgncid, inh,
                                                          candidate_variants, inheritance_report)
                elif inh == 'Mosaic':
                    mosaic_heterozygous_parents_filter(v,
                                                          variants[v]['child'],
                                                          hgncid, inh,
                                                          candidate_variants, inheritance_report)
                elif inh == 'Imprinted':
                    imprinted_heterozygous_parents_filter(v, variants[v]['child'],
                                                     hgncid, inh, candidate_variants, inheritance_report)
                else:
                    logging.info(v + " unknown gene mode " + inh)
        elif variants[v]['child'].genotype == '2':
            #homozygous
            for inh in gene['mode']:
                if inh == 'Biallelic':
                    biallelic_homozygous_parents_filter(v,
                                                          variants[v]['child'],
                                                          hgncid, inh,
                                                          candidate_variants, inheritance_report)
                elif inh == 'Monoallelic':
                    monoallelic_homozygous_parents_filter(v,
                                                            variants[v][
                                                                'child'],
                                                            hgncid, inh,
                                                            candidate_variants, inheritance_report)
                elif inh == 'Mosaic':
                    mosaic_homozygous_parents_filter(v,
                                                       variants[v]['child'],
                                                       hgncid, inh,
                                                       candidate_variants, inheritance_report)
                elif inh == 'Imprinted':
                    imprinted_homozygous_parents_filter(v,
                                                          variants[v]['child'],
                                                          hgncid, inh,
                                                          candidate_variants, inheritance_report)
                else:
                    logging.info(v + " unknown gene mode " + inh)
        else:
            logging.info(v + " fails inheritance filters - child must be hom or"
                             " het")



def autosomal_no_parents(hgncid, gene, variants, family, candidate_variants, inheritance_report,
                           lof_cqs):
    '''screens variants in autosomes where there are both parents and adds
    candidates to candidate_variants
    args:
    gene - info about the gene
    variants - all variants in this gene which pass preinheritance filters
    candidate_variants - output dict
    '''
    pass

def autosomal_single_parent(hgncid, gene, variants, family, candidate_variants, inheritance_report,
                            lof_cqs):
    '''screens variants in autosomes where there is one parent and adds
    candidates to candidate_variants
    args:
    gene - info about the gene
    variants - all variants in this gene which pass preinheritance filters
    candidate_variants - output dict
    '''
    pass

def biallelic_heterozygous_parents_filter(varid, var, hgncid, inh, candidate_variants, inheritance_report):
    pass

def monoallelic_heterozygous_parents_filter(varid, var, hgncid, inh, candidate_variants, inheritance_report):
    pass

def mosaic_heterozygous_parents_filter(varid, var, hgncid, inh, candidate_variants, inheritance_report):
    pass

def imprinted_heterozygous_parents_filter(varid, var, hgncid, inh, candidate_variants, inheritance_report):
    pass

def biallelic_homozygous_parents_filter(varid, var, hgncid, inh, candidate_variants, inheritance_report):
    print('here')
    exit(0)

def monoallelic_homozygous_parents_filter(varid, var, hgncid, inh, candidate_variants, inheritance_report):
    pass

def mosaic_homozygous_parents_filter(varid, var, hgncid, inh, candidate_variants, inheritance_report):
    pass

def imprinted_homozygous_parents_filter(varid, var, hgncid, inh, candidate_variants, inheritance_report):
    pass

