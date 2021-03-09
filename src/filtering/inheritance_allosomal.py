import logging

from utils.utils import common_elements
from utils.utils import add_compound_het_to_candidates
from utils.utils import add_single_var_to_candidates


def allosomal_both_parents(hgncid, gene, variants, family, candidate_variants,
                           lof_cqs):
    '''screens variants in X/Y where there are both parents and adds candidates
    to candidate_variants
    args:
    gene - info about the gene
    variants - all variants in this gene which pass preinheritance filters
    candidate_variants - output dict
    '''
    for v in variants.keys():
        cqs = variants[v]['child'].consequence.split('&')
        cqs_in_lof = common_elements(cqs, lof_cqs)
        is_lof = False
        if len(cqs_in_lof) > 0:
            is_lof = True
        if variants[v]['child'].genotype == "1":
            # heterozygous
            for inh in gene['mode']:
                if inh == 'X-linked dominant' or \
                        inh == 'X-linked over-dominance' or inh == 'Hemizygous':
                    if variants[v]['child'].is_dad_hom_ref() and \
                            variants[v]['child'].is_mum_hom_ref():
                        add_single_var_to_candidates(v, variants[v]['child'],
                                                     hgncid, inh,
                                                     candidate_variants)

                    elif (family.mum.get_affected_status() and \
                          not variants[v]['child'].is_mum_hom_ref()) and \
                            (family.dad.get_affected_status() or \
                             variants[v]['child'].is_dad_hom_ref()):
                        # inherited from affected mum, dad either affected or hom ref
                        if inh == 'Hemizygous':
                            add_compound_het_to_candidates(v,
                                                           variants[v]['child'],
                                                           hgncid, inh,
                                                           candidate_variants)

                        else:
                            add_single_var_to_candidates(v,
                                                         variants[v]['child'],
                                                         hgncid, inh,
                                                         candidate_variants)

                    elif (family.dad.get_affected_status() and \
                          not variants[v]['child'].is_dad_hom_ref()) and \
                            (family.mum.get_affected_status() or \
                             variants[v]['child'].is_mum_hom_ref()):
                        # inherited from affected dad, mum either affected or hom ref
                        if inh == 'Hemizygous':
                            add_compound_het_to_candidates(v,
                                                           variants[v]['child'],
                                                           hgncid, inh,
                                                           candidate_variants)

                        else:
                            add_single_var_to_candidates(v,
                                                         variants[v]['child'],
                                                         hgncid, inh,
                                                         candidate_variants)

                    elif is_lof and inh == 'X-linked over-dominance' and \
                            not variants[v]['child'].is_dad_hom_ref():
                        add_single_var_to_candidates(v, variants[v]['child'],
                                                     hgncid, inh,
                                                     candidate_variants)

                    else:
                        logging.info(v + " failed allosomal inheritance filter "
                                         "heterozgous variant")
                # elif inh == 'Hemizygous':
                #     pass
                else:
                    logging.info(v + " failed allosomal inheritance filter - "
                                     "invalid inheritance type " + inh)
        elif variants[v]['child'].genotype == "2":
            # homozygous
            if family.proband.gender == 'M':
                # Male proband
                for inh in gene['mode']:
                    if inh in ['X-linked dominant', 'Hemizygous',
                               'X-linked over-dominance']:
                        if variants[v]['child'].is_mum_hom_ref():
                            # mum hom ref
                            add_single_var_to_candidates(v,
                                                         variants[v]['child'],
                                                         hgncid, inh,
                                                         candidate_variants)

                        elif inh == 'X-linked over-dominance' and \
                                variants[v]['child'].is_mum_het() and \
                                not family.mum.get_affected_status():
                            logging.info(v + " failed allosomal inheritance "
                                             "filter hom in male inherited from"
                                             " unaffected mother in X-linked "
                                             "over-dominant gene")
                        elif variants[v]['child'].is_mum_het() or \
                                (variants[v]['child'].is_mum_hom_alt() and \
                                 family.mum.get_affected_status()):
                            add_single_var_to_candidates(v,
                                                         variants[v]['child'],
                                                         hgncid, inh,
                                                         candidate_variants)

                        else:
                            logging.info(v + " failed allosomal inheritance "
                                             "filter homozygous variant in male"
                                             "proband")
                    else:
                        logging.info(v + " failed allosomal inheritance filter "
                                         "- invalid inheritance type " + inh)
            elif family.proband.gender == 'F':
                # Female proband
                for inh in gene['mode']:
                    if inh in ['X-linked dominant', 'Hemizygous',
                               'X-linked over-dominance']:
                        if variants[v]['child'].is_mum_hom_ref() or \
                                variants[v]['child'].is_dad_hom_ref():
                            logging.info(
                                v + " failed allosomal inheritance filter "
                                    "- homozygous variant in female "
                                    "proband and parent(s) hom ref")
                        elif (variants[v]['child'].is_mum_het() or \
                              (variants[v]['child'].is_mum_hom_alt() and \
                               family.mum.get_affected_status())) and \
                                (variants[v]['child'].is_dad_hom_alt() and \
                                 family.dad.get_affected_status()):
                            add_single_var_to_candidates(v,
                                                         variants[v]['child'],
                                                         hgncid, inh,
                                                         candidate_variants)

                        else:
                            logging.info(
                                v + " failed allosomal inheritance filter "
                                    "- homozygous variant in female")
                    else:
                        logging.info(v + " failed allosomal inheritance filter "
                                         "- invalid inheritance type " + inh)

            else:
                logging.info(v + " failed allosomal inheritance filter unknown "
                                 "proband gender " + family.proband.gender)
        else:
            logging.info(v + " failed allosomal inheritance filter, unknown "
                             "genotype")


def allosomal_no_parents(hgncid, gene, variants, candidate_variants, lof_cqs):
    '''screens variants in X/Y where there are no parents and adds candidates
    to candidate_variants
    args:
    gene - info about the gene
    variants - all variants in this gene which pass preinheritance filters
    candidate_variants - output dict
    '''
    for v in variants.keys():
        for inh in gene['mode']:
            if inh == 'X-linked dominant' or \
                    inh == 'X-linked over-dominance' or \
                    inh == 'Hemizygous':
                add_single_var_to_candidates(v, variants[v]['child'], hgncid,
                                             inh, candidate_variants)
            else:
                logging.info(v + " failed allosomal inheritance filter no "
                                 "parents")


def allosomal_single_parent(hgncid, gene, variants_per_gene, family,
                            candidate_variants, lof_cqs):
    '''screens variants in X/Y where there is one parent and adds candidates
    to candidate_variants
    args:
    gene - info about the gene
    variants - all variants in this gene which pass preinheritance filters
    candidate_variants - output dict
    '''
    pass
