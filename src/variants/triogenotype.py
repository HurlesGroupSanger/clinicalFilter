"""
copyright
"""
import logging
from variants.variant import Variant

def add_trio_genotypes(family, variants):
    '''add trio genotypes to all child variants for a family'''
    if family.has_both_parents():
        #if parents are present we assume positions not present are ref/ref
        add_trio_genotypes_both_parents(variants)
    elif family.has_mum():
        add_trio_genotypes_mum_only(variants)
    elif family.has_dad():
        add_trio_genotypes_dad_only(variants)
    elif family.has_no_parents():
        add_trio_genotypes_no_parents(variants)
    else:
        print(
            "Error: should not get here as family must either have both, one \
             or no parents")
        logging.error("Can't calculate trio genotypes - family error")
        exit(1)

def add_trio_genotypes_both_parents(variants):
    '''add trio genotypes in a dict of variants where there are both parents'''
    for v in variants['child'].keys():
        childgeno = variants['child'][v].genotype
        mumgeno = '0'
        dadgeno = '0'
        if v in variants['mum'].keys():
            mumgeno = variants['mum'][v].genotype
        if v in variants['dad'].keys():
            dadgeno = variants['dad'][v].genotype

        triogenotype = childgeno + mumgeno + dadgeno
        variants['child'][v].set_triogenotype(triogenotype)

def add_trio_genotypes_no_parents(variants):
    '''add trio genotypes in a dict of variants where there are no parents'''
    for v in variants['child'].keys():
        childgeno = variants['child'][v].genotype
        mumgeno = 'NA'
        dadgeno = 'NA'
        triogenotype = childgeno + mumgeno + dadgeno
        variants['child'][v].set_triogenotype(triogenotype)

def add_trio_genotypes_mum_only(variants):
    '''add trio genotypes in a dict of variants where there is mum only'''
    for v in variants['child'].keys():
        childgeno = variants['child'][v].genotype
        mumgeno = '0'
        dadgeno = 'NA'
        if v in variants['mum'].keys():
            mumgeno = variants['mum'][v].genotype

        triogenotype = childgeno + mumgeno + dadgeno
        variants['child'][v].set_triogenotype(triogenotype)

def add_trio_genotypes_dad_only(variants):
    '''add trio genotype in a dict of variants where there is dad only'''
    for v in variants['child'].keys():
        childgeno = variants['child'][v].genotype
        mumgeno = 'NA'
        dadgeno = '0'
        if v in variants['dad'].keys():
            dadgeno = variants['dad'][v].genotype

        triogenotype = childgeno + mumgeno + dadgeno
        variants['child'][v].set_triogenotype(triogenotype)



