"""copyright"""

import logging
from itertools import combinations

from filtering.inheritance_autosomal import autosomal_no_parents
from filtering.inheritance_autosomal import autosomal_both_parents
from filtering.inheritance_autosomal import autosomal_single_parent
from filtering.inheritance_allosomal import allosomal_no_parents
from filtering.inheritance_allosomal import allosomal_both_parents
from filtering.inheritance_allosomal import allosomal_single_parent


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
        for hgncid in variants_per_gene.keys():
            if hgncid in genes.keys():
                if genes[hgncid]['chr'] in ['X', 'Y']:
                    if parents == 'both':
                        allosomal_both_parents(hgncid, genes[hgncid],
                                               variants_per_gene[hgncid],
                                               family, candidate_variants,
                                               lof_cqs)
                    elif parents == 'none':
                        allosomal_no_parents(hgncid, genes[hgncid],
                                             variants_per_gene[hgncid],
                                             candidate_variants, lof_cqs)
                    else:
                        allosomal_single_parent(hgncid, genes[hgncid],
                                                variants_per_gene[hgncid],
                                                family, candidate_variants,
                                                lof_cqs)
                else:
                    if parents == 'both':
                        autosomal_both_parents(hgncid, genes[hgncid],
                                               variants_per_gene[hgncid],
                                               family, candidate_variants,
                                               lof_cqs)
                    elif parents == 'none':
                        autosomal_no_parents(hgncid, genes[hgncid],
                                             variants_per_gene[hgncid],
                                             candidate_variants, lof_cqs)
                    else:
                        autosomal_single_parent(hgncid, genes[hgncid],
                                                variants_per_gene[hgncid],
                                                family, candidate_variants,
                                                lof_cqs)
            else:
                for v in variants_per_gene[hgncid].keys():
                    logging.info(
                        v + " gene not in DDG2P: " +
                        variants_per_gene[hgncid][v]['child'].symbol)

    if regions:
        # todo filtering for regions
        pass

    if trusted_variants:
        # todo filtering for trusted variants
        pass

    screen_compound_hets(candidate_variants, family)
    return candidate_variants


def screen_compound_hets(candidate_variants, family):
    # sort out compound hets
    compound_het_passes = {}
    for gn in candidate_variants['compound_hets'].keys():
        # print(gn)
        if len(candidate_variants['compound_hets'][gn].keys()) < 2:
            for v in candidate_variants['compound_hets'][gn].keys():
                logging.info(
                    v + " failed compound het screen: <2 vars in hgnc " + gn)

        combs_to_screen = list(
            combinations(candidate_variants['compound_hets'][gn].keys(), 2))
        for pair in combs_to_screen:
            var1 = candidate_variants['compound_hets'][gn][pair[0]]['variant']
            var2 = candidate_variants['compound_hets'][gn][pair[1]]['variant']
            if var1 != var2:
                if is_compound_pair(pair[0], var1, pair[1], var2, family):
                    if not gn in compound_het_passes.keys():
                        compound_het_passes[gn] = {}
                    if not pair[0] in compound_het_passes[gn].keys():
                        compound_het_passes[gn][pair[0]] = {}
                    if not pair[1] in compound_het_passes[gn].keys():
                        compound_het_passes[gn][pair[1]] = {}

                    compound_het_passes[gn][pair[0]]['variant'] = var1
                    compound_het_passes[gn][pair[0]]['mode'] = \
                    candidate_variants['compound_hets'][gn][pair[0]]['mode']
                    compound_het_passes[gn][pair[1]]['variant'] = var2
                    compound_het_passes[gn][pair[1]]['mode'] = \
                        candidate_variants['compound_hets'][gn][pair[1]]['mode']

    candidate_variants['compound_hets'] = compound_het_passes
    # for gn in candidate_variants['compound_hets'].keys():
    #     print(gn)
    # print(candidate_variants['compound_hets'])
    pass


def is_compound_pair(varid1, var1, varid2, var2, family):
    '''Test to see if a pair of variants could be a compound het'''

    if family.has_no_parents():
        """If there are no parents and the variants are not both missense then 
        pass"""
        if var1.consequence == 'missense_variant' and \
                var2.consequence == 'missense_variant':
            logging.info(varid1 + " " + varid2 + " failed compound het "
                                                 "screen, no parents and both missense")
            return False
        elif (var1.pid != '.' and var2.pid != '.') and (var1.pid == var2.pid):
            logging.info(varid1 + " " + varid2 + " failed compound het "
                                                 "screen, variants in cis")
            return False

        else:
            return True
    elif family.has_both_parents():
        if var1.chrom == 'X' and not family.dad.get_affected_status() and \
                (var1.get_dad_genotype() == '0' or \
                 var2.get_dad_genotype() == '0'):
            logging.info(varid1 + " " + varid2 + " failed compound het "
                                                 "screen, X chrom and dad unaffected and hom ref for"
                                                 " 1 variant")
            return False
        elif (var1.get_mum_genotype() == '0' and \
              var2.get_mum_genotype() != '0' and \
              var1.get_dad_genotype() != '0' and \
              var2.get_dad_genotype() == '0') or \
                (var2.get_mum_genotype() == '0' and \
                 var1.get_mum_genotype() != '0' and \
                 var2.get_dad_genotype() != '0' and \
                 var1.get_dad_genotype() == '0'):
            '''one variant is inherited from each parent'''
            return True
        elif ((var1.denovo_snv or var1.denovo_indel) and \
              var2.get_mum_genotype() != '0' and \
              var2.get_dad_genotype() == '0') or \
                ((var1.denovo_snv or var1.denovo_indel) and \
                 var2.get_mum_genotype() == '0' and \
                 var2.get_dad_genotype() != '0') or \
                ((var2.denovo_snv or var2.denovo_indel) and \
                 var1.get_mum_genotype() != '0' and \
                 var1.get_dad_genotype() == '0') or \
                ((var2.denovo_snv or var2.denovo_indel) and \
                 var1.get_mum_genotype() == '0' and \
                 var1.get_dad_genotype() != '0'):
            '''one variant is DNM and the other is inherited'''
            return True
        else:
            logging.info(varid1 + " " + varid2 + " failed compound het screen")
            return False

    elif family.has_dad():
        pass
    # todo screening of compound hets with one parent
    elif family.has_mum():
        pass
    else:
        print(
            "Error: should not get here as family must either have both, one \
             or no parents")
        logging.error("Can't parse compound hets - family error")
        exit(1)
    logging.info(varid1 + " " + varid2 + " failed compound het screen")
    return False
