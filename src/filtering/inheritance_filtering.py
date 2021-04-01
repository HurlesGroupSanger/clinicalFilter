"""copyright"""

import logging

from filtering.inheritance_autosomal import autosomal_no_parents
from filtering.inheritance_autosomal import autosomal_both_parents
from filtering.inheritance_autosomal import autosomal_single_parent
from filtering.inheritance_allosomal import allosomal_no_parents
from filtering.inheritance_allosomal import allosomal_both_parents
from filtering.inheritance_allosomal import allosomal_single_parent
from filtering.compound_hets import screen_compound_hets


def inheritance_filter(variants_per_gene, family, genes, regions,
                       trusted_variants):
    # Autosomal without parents
    # Autosomal with parents
    #   Heterozygous
    #   Homozygous
    # Autosomal with one parent
    # Allosomal without parents
    # Allosomal with parents
    #   Heterozygous
    #   Homozygous
    # Allosomal with one parent
    # Identification of compound hets,

    candidate_variants = {}
    candidate_variants['single_variants'] = {}
    candidate_variants['compound_hets'] = {}

    # inheritance report gives data to populate a matrix for every possible
    # combination of proband and parent GT and affected status for each mode of
    # inheritance. Some variants will be counted twice (genes which are both
    # mono and biallelic). Split into X and autosome
    inheritance_report = {}

    lof_cqs = ['transcript_ablation', 'splice_donor_variant', 'stop_lost',
               'splice_acceptor_variant', 'stop_gained', 'frameshift_variant',
               'start_lost']

    parents = ''
    if family.has_both_parents():
        parents = 'both'
    elif family.has_no_parents():
        parents = 'none'
    elif family.has_dad():
        parents = 'dad_only'
    elif family.has_mum():
        parents = 'mum_only'
    else:
        print(
            "Error: should not get here as family must either have both, one \
             or no parents")
        logging.error("Can't calculate inheritance type - family error")
        exit(1)

    if genes:
        inheritance_filter_genes(variants_per_gene, family, genes, parents,
                                 candidate_variants, inheritance_report,
                                 lof_cqs)

    if regions:
        # todo filtering for regions
        pass

    if trusted_variants:
        # todo filtering for trusted variants
        pass

    screen_compound_hets(candidate_variants, family)
    return candidate_variants, inheritance_report

def inheritance_filter_genes(variants_per_gene, family, genes, parents, candidate_variants, inheritance_report, lof_cqs):
    # inheritance filters for use with a gene list
    for hgncid in variants_per_gene.keys():
        if hgncid in genes.keys():
            if genes[hgncid]['chr'] in ['X', 'Y']:
                if parents == 'both':
                    allosomal_both_parents(hgncid, genes[hgncid],
                                           variants_per_gene[hgncid],
                                           family, candidate_variants,
                                           inheritance_report,
                                           lof_cqs)
                elif parents == 'none':
                    allosomal_no_parents(hgncid, genes[hgncid],
                                         variants_per_gene[hgncid],
                                         candidate_variants, inheritance_report,
                                         lof_cqs)
                else:
                    allosomal_single_parent(hgncid, genes[hgncid],
                                            variants_per_gene[hgncid],
                                            family, candidate_variants,
                                            inheritance_report,
                                            lof_cqs)
            else:
                if parents == 'both':
                    autosomal_both_parents(hgncid, genes[hgncid],
                                           variants_per_gene[hgncid],
                                           family, candidate_variants,
                                           inheritance_report,
                                           lof_cqs)
                elif parents == 'none':
                    autosomal_no_parents(hgncid, genes[hgncid],
                                         variants_per_gene[hgncid],
                                         candidate_variants, inheritance_report,
                                         lof_cqs)
                else:
                    autosomal_single_parent(hgncid, genes[hgncid],
                                            variants_per_gene[hgncid],
                                            family, candidate_variants,
                                            inheritance_report,
                                            lof_cqs)
        else:
            for v in variants_per_gene[hgncid].keys():
                logging.info(
                    v + " gene not in DDG2P: " +
                    variants_per_gene[hgncid][v]['child'].symbol)
