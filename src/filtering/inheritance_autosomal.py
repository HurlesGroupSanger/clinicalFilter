import logging

from utils.utils import common_elements
from utils.utils import add_compound_het_to_candidates
from utils.utils import add_single_var_to_candidates

def autosomal_no_parents(hgncid, gene, variants, candidate_variants, lof_cqs):
    '''screens variants in autosomes where there are no parents and adds
    candidates to candidate_variants
    args:
    gene - info about the gene
    variants - all variants in this gene which pass preinheritance filters
    candidate_variants - output dict
    '''

    for v in variants.keys():
        # is there a lof cq?
        cqs = variants[v]['child'].consequence.split('&')
        cqs_in_lof = common_elements(cqs, lof_cqs)
        is_lof = False
        if len(cqs_in_lof) > 0:
            is_lof = True

        for inh in gene['mode']:
            if inh == 'Biallelic':
                if variants[v]['child'].genotype == "1":
                    add_compound_het_to_candidates(v, variants[v]['child'], hgncid, inh, candidate_variants)
                elif variants[v]['child'].genotype == "2":
                    add_single_var_to_candidates(v, variants[v]['child'], hgncid, inh, candidate_variants)
                else:
                    logging.info(v + " fails inheritance filters - child must "
                                     "be hom or het")

            elif inh == 'Monoallelic' and variants[v]['child'].genotype == "1":
                add_single_var_to_candidates(v, variants[v]['child'], hgncid, inh,
                                             candidate_variants)

            elif inh == 'Imprinted' and variants[v]['child'].genotype == "1":
                if is_lof:
                    add_single_var_to_candidates(v, variants[v]['child'], hgncid, inh,
                                                 candidate_variants)

                else:
                    logging.info(v + " fails inheritance filters for imprinted "
                                     "- child must be LoF")
            else:
                logging.info(v + "failed heterozygous inheritance filters "
                                 "invalid inheritance type" + inh)


def autosomal_both_parents(hgncid, gene, variants, family, candidate_variants,
                           lof_cqs):
    '''screens variants in autosomes where there are both parents and adds
    candidates to candidate_variants
    args:
    gene - info about the gene
    variants - all variants in this gene which pass preinheritance filters
    candidate_variants - output dict
    '''
    for v in variants.keys():
        # is there a lof cq?
        cqs = variants[v]['child'].consequence.split('&')
        cqs_in_lof = common_elements(cqs, lof_cqs)
        is_lof = False
        if len(cqs_in_lof) > 0:
            is_lof = True
        if variants[v]['child'].genotype == "1":
            # heterozygous
            rejected_variants = {}
            not_single_var = {}
            parents_hom_ref = False
            transmitted_from_affected_parent = False
            if not variants[v]['child'].triogenotype == "100":
                parents_hom_ref = True
            if family.dad.get_affected_status() and \
                    not variants[v]['child'].is_dad_hom_ref() and \
                    (family.mum.get_affected_status() or
                     variants[v]['child'].is_dad_hom_ref()):
                transmitted_from_affected_parent = True
            if family.mum.get_affected_status() and \
                    not variants[v]['child'].is_mum_hom_ref() and \
                    (family.dad.get_affected_status() or
                     variants[v]['child'].is_dad_hom_ref()):
                transmitted_from_affected_parent = True

            for inh in gene['mode']:
                # consider each inheritance mode separately
                if inh == 'Monoallelic' or inh == 'Mosaic':
                    if parents_hom_ref or transmitted_from_affected_parent:
                        add_single_var_to_candidates(v, variants[v]['child'], hgncid, inh,
                                                     candidate_variants)
                    else:
                        logging.info(v + "failed heterozygous inheritance "
                                         "filters " + inh + " and parents not "
                                                            "hom ref or not transmitted from "
                                                            "affected parent")
                elif inh == 'Biallelic':
                    if not parents_hom_ref and not transmitted_from_affected_parent:
                        if (not variants[v]['child'].is_dad_hom_alt() or \
                            family.dad.get_affected_status() == True) and \
                                (not variants[v]['child'].is_mum_hom_alt() or \
                                 family.mum.get_affected_status() == True):
                            add_compound_het_to_candidates(v, variants[v]['child'],
                                                           hgncid, inh,
                                                           candidate_variants)

                        else:
                            logging.info(v + "failed heterozygous inheritance "
                                             "filters - inherited from hom alt "
                                             "unaffcted parent")
                    else:
                        add_compound_het_to_candidates(v, variants[v]['child'],
                                                       hgncid, inh,
                                                       candidate_variants)

                elif inh == 'Imprinted':
                    if is_lof:
                        if parents_hom_ref or transmitted_from_affected_parent:
                            add_single_var_to_candidates(v, variants[v]['child'], hgncid, inh,
                                                         candidate_variants)
                        else:
                            logging.info(v + "failed heterozygous inheritance "
                                             "filters, Imprinted var and parents "
                                             "not hom ref or not transmitted from "
                                             "affected parent")
                    else:
                        logging.info(v + "failed heterozygous inheritance "
                                         " Imprinted variants must be LoF")

                else:
                    logging.info(v + "failed heterozygous inheritance filters -"
                                     " invalid inheritance type" + inh)


        elif variants[v]['child'].genotype == "2":
            # homozygous
            if variants[v]['child'].is_mum_hom_ref() or \
                    variants[v]['child'].is_dad_hom_ref():
                logging.info(v + " homozygous variant failed due to hom ref "
                                 "parent(s)")
            #     todo CNVs will change this
            else:
                for inh in gene['mode']:
                    # consider each inheritance mode separately
                    if inh == 'Biallelic':
                        if variants[v]['child'].is_mum_het() and \
                                variants[v]['child'].is_dad_het():
                            add_single_var_to_candidates(v, variants[v]['child'], hgncid, inh,
                                                         candidate_variants)

                        elif (variants[v]['child'].is_mum_hom_alt() and \
                              family.mum.get_affected_status()) and \
                                (not variants[v]['child'].is_dad_hom_alt() or \
                                 family.dad.get_affected_status()):
                            add_single_var_to_candidates(v, variants[v]['child'], hgncid, inh,
                                                         candidate_variants)

                        elif (variants[v]['child'].is_dad_hom_alt() and \
                              family.dad.get_affected_status()) and \
                                (not variants[v]['child'].is_mum_hom_alt() or \
                                 family.mum.get_affected_status()):
                            add_single_var_to_candidates(v, variants[v]['child'], hgncid, inh,
                                                         candidate_variants)

                        else:
                            logging.info(
                                v + " homozygous biallelic variant failed "
                                    " homozygous parents not affected")

                    elif inh == 'Monoallelic':
                        if family.mum.get_affected_status() and \
                                family.dad.get_affected_status():
                            add_single_var_to_candidates(v, variants[v]['child'], hgncid, inh,
                                                         candidate_variants)

                        else:
                            logging.info(v + " homozygous monoallelic variant "
                                             "failed - not inherited from "
                                             "affected parents")
                    elif inh == 'Imprinted':
                        if is_lof:
                            if not variants[v]['child'].is_dad_hom_ref() or \
                                not variants[v]['child'].is_mum_hom_ref():
                                add_single_var_to_candidates(v, variants[v]['child'],
                                                             hgncid, inh,
                                                             candidate_variants)

                            else:
                                logging.info(v + "failed homozygous inheritance"
                                                 " imprinted variants must have "
                                                 "hom ref parents")
                        else:
                            logging.info(v + "failed homozygous inheritance "
                                             " Imprinted variants must be LoF")
                    else:
                        logging.info(
                            v + "failed homozygous inheritance filters -"
                                " invalid inheritance type" + inh)


        else:
            print(
                "Error: should not get here as genotype must be het (1) or "
                "hom (2)")
            logging.error("Can't calculate inheritance type - family error")
            exit(1)
        # now check if the variant is in candidates, and if not print to log file

        if v not in candidate_variants['single_variants'].keys():
            # for any variants not in candidate variants print to log file
            if hgncid not in candidate_variants['compound_hets'].keys():
                logging.info(v + " failed autosomal inheritance filters (both "
                                 "parents present) for " + gene['symbol'])
            else:
                if v not in candidate_variants['compound_hets'][hgncid].keys():
                    logging.info(
                        v + " failed autosomal inheritance filters (both "
                            "parents present) for " + gene['symbol'])


def autosomal_single_parent(hgncid, gene, variants, family, candidate_variants,
                            lof_cqs):
    '''screens variants in autosomes where there is one parent and adds
    candidates to candidate_variants
    args:
    gene - info about the gene
    variants - all variants in this gene which pass preinheritance filters
    candidate_variants - output dict
    '''
    pass
