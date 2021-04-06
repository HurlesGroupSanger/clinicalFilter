import logging

from utils.utils import common_elements
from utils.utils import add_compound_het_to_candidates
from utils.utils import add_single_var_to_candidates
from utils.utils import convert_genotype_to_gt
from filtering.inheritance_report import populate_inheritance_report

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
                    biallelic_heterozygous_parents_filter(v, variants[v]['child'], family,
                                                     hgncid, candidate_variants, inheritance_report)
                elif inh == 'Monoallelic':
                    monoallelic_heterozygous_parents_filter(v,
                                                          variants[v]['child'], family,
                                                          hgncid,
                                                          candidate_variants, inheritance_report)
                elif inh == 'Mosaic':
                    mosaic_heterozygous_parents_filter(v,
                                                          variants[v]['child'], family,
                                                          hgncid,
                                                          candidate_variants, inheritance_report)
                elif inh == 'Imprinted':
                    imprinted_heterozygous_parents_filter(v, variants[v]['child'], family,
                                                     hgncid, candidate_variants, inheritance_report)
                else:
                    logging.info(v + " unknown gene mode " + inh)
        elif variants[v]['child'].genotype == '2':
            #homozygous
            for inh in gene['mode']:
                if inh == 'Biallelic':
                    biallelic_homozygous_parents_filter(v,
                                                          variants[v]['child'], family,
                                                          hgncid,
                                                          candidate_variants, inheritance_report)
                elif inh == 'Monoallelic':
                    monoallelic_homozygous_parents_filter(v,
                                                            variants[v][
                                                                'child'], family,
                                                            hgncid,
                                                            candidate_variants, inheritance_report)
                elif inh == 'Mosaic':
                    mosaic_homozygous_parents_filter(v,
                                                       variants[v]['child'], family,
                                                       hgncid,
                                                       candidate_variants, inheritance_report)
                elif inh == 'Imprinted':
                    imprinted_homozygous_parents_filter(v,
                                                          variants[v]['child'], family,
                                                          hgncid,
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

def biallelic_heterozygous_parents_filter(varid, var, family, hgncid, candidate_variants, inheritance_report):
    child_gt = var.gt
    mum_genotype = var.get_mum_genotype()
    dad_genotype = var.get_mum_genotype()
    mum_gt = convert_genotype_to_gt(mum_genotype)
    dad_gt = convert_genotype_to_gt(dad_genotype)
    mum_aff = family.mum.affected
    dad_aff = family.dad.affected

    populate_inheritance_report(inheritance_report, 'autosomal', 'biallelic', child_gt, mum_gt, dad_gt, mum_aff, dad_aff)
    vpass = 'n'
    if mum_aff and dad_aff:
        if not dad_gt == '1/1' and not mum_gt == '1/1':
            add_compound_het_to_candidates(varid, var, hgncid, 'biallelic',
                                       candidate_variants)
            vpass = 'y'
    elif mum_aff and not dad_aff:
        if not dad_gt == '1/1':
            add_compound_het_to_candidates(varid, var, hgncid, 'biallelic',
                                           candidate_variants)
            vpass = 'y'
    elif not mum_aff and dad_aff:
        if not mum_gt == '1/1':
            add_compound_het_to_candidates(varid, var, hgncid, 'biallelic',
                                           candidate_variants)
            vpass = 'y'
    else:
        if not dad_gt == '1/1' and not mum_gt == '1/1':
            add_compound_het_to_candidates(varid, var, hgncid, 'biallelic',
                                           candidate_variants)
            vpass = 'y'

    if vpass == 'n':
        logging.info(varid + " failed inheritance filter for heterozygous "
                             "variant in biallelic gene")

def monoallelic_heterozygous_parents_filter(varid, var, family,  hgncid, candidate_variants, inheritance_report):

    child_gt = var.gt
    mum_genotype = var.get_mum_genotype()
    dad_genotype = var.get_mum_genotype()
    mum_gt = convert_genotype_to_gt(mum_genotype)
    dad_gt = convert_genotype_to_gt(dad_genotype)
    mum_aff = family.mum.affected
    dad_aff = family.dad.affected

    populate_inheritance_report(inheritance_report, 'autosomal', 'monoallelic',
                                child_gt, mum_gt, dad_gt, mum_aff, dad_aff)
    vpass = 'n'

    if mum_aff and dad_aff:
        if not dad_gt == '1/1' and not mum_gt == '1/1':
            add_single_var_to_candidates(varid, var, hgncid, 'monoallelic', candidate_variants)
            vpass = 'y'
    elif mum_aff and not dad_aff:
        if dad_gt == '0/0':
            add_single_var_to_candidates(varid, var, hgncid, 'monoallelic', candidate_variants)
            vpass = 'y'
    elif not mum_aff and dad_aff:
        if mum_gt == '0/0':
            add_single_var_to_candidates(varid, var, hgncid, 'monoallelic', candidate_variants)
            vpass = 'y'
    else:
        if mum_gt == '0/0' and dad_gt == '0/0':
            add_single_var_to_candidates(varid, var, hgncid, 'monoallelic', candidate_variants)
            vpass = 'y'

    if vpass == 'n':
        logging.info(varid + " failed inheritance filter for heterozygous "
                             "variant in monoallelic gene")



def mosaic_heterozygous_parents_filter(varid, var, family, hgncid, candidate_variants, inheritance_report):

    child_gt = var.gt
    mum_genotype = var.get_mum_genotype()
    dad_genotype = var.get_mum_genotype()
    mum_gt = convert_genotype_to_gt(mum_genotype)
    dad_gt = convert_genotype_to_gt(dad_genotype)
    mum_aff = family.mum.affected
    dad_aff = family.dad.affected

    populate_inheritance_report(inheritance_report, 'autosomal', 'mosaic',
                                child_gt, mum_gt, dad_gt, mum_aff, dad_aff)
    vpass = 'n'

    if mum_aff and dad_aff:
        if not dad_gt == '1/1' and not mum_gt == '1/1':
            add_single_var_to_candidates(varid, var, hgncid, 'mosaic', candidate_variants)
            vpass = 'y'
    elif mum_aff and not dad_aff:
        if dad_gt == '0/0':
            add_single_var_to_candidates(varid, var, hgncid, 'mosaic', candidate_variants)
            vpass = 'y'
    elif not mum_aff and dad_aff:
        if mum_gt == '0/0':
            add_single_var_to_candidates(varid, var, hgncid, 'mosaic', candidate_variants)
            vpass = 'y'
    else:
        if mum_gt == '0/0' and dad_gt == '0/0':
            add_single_var_to_candidates(varid, var, hgncid, 'mosaic', candidate_variants)
            vpass = 'y'

    if vpass == 'n':
        logging.info(varid + " failed inheritance filter for heterozygous "
                             "variant in mosaic gene")

def imprinted_heterozygous_parents_filter(varid, var, family, hgncid, candidate_variants, inheritance_report):
    pass

def biallelic_homozygous_parents_filter(varid, var, family, hgncid, candidate_variants, inheritance_report):
    # todo will need modification when CNVs (and UPDs) added
    child_gt = var.gt
    print(child_gt)
    mum_genotype = var.get_mum_genotype()
    dad_genotype = var.get_mum_genotype()
    mum_gt = convert_genotype_to_gt(mum_genotype)
    dad_gt = convert_genotype_to_gt(dad_genotype)
    mum_aff = family.mum.affected
    dad_aff = family.dad.affected

    populate_inheritance_report(inheritance_report, 'autosomal', 'biallelic',
                                child_gt, mum_gt, dad_gt, mum_aff, dad_aff)
    vpass = 'n'

    if mum_aff and dad_aff:
        if mum_gt != '0/0' and dad_gt != '0/0':
            add_single_var_to_candidates(varid, var, hgncid, 'biallelic',
                                         candidate_variants)
            vpass = 'y'
    elif mum_aff and not dad_aff:
        if dad_gt == '0/1' and mum_gt != '0/0':
            add_single_var_to_candidates(varid, var, hgncid, 'biallelic',
                                         candidate_variants)
            vpass = 'y'
    elif not mum_aff and dad_aff:
        if mum_gt == '0/1' and dad_gt != '0/0':
            add_single_var_to_candidates(varid, var, hgncid, 'biallelic',
                                         candidate_variants)
            vpass = 'y'
    else:
        if mum_gt == '0/1' and dad_gt == '0/1':
            add_single_var_to_candidates(varid, var, hgncid, 'biallelic',
                                         candidate_variants)
            vpass = 'y'

    if vpass == 'n':
        logging.info(varid + " failed inheritance filter for heterozygous "
                             "variant in mosaic gene")

def monoallelic_homozygous_parents_filter(varid, var, family, hgncid, candidate_variants, inheritance_report):
    pass

def mosaic_homozygous_parents_filter(varid, var, family, hgncid, candidate_variants, inheritance_report):
    pass

def imprinted_homozygous_parents_filter(varid, var, family, hgncid, candidate_variants, inheritance_report):
    pass

