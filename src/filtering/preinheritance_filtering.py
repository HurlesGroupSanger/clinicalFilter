"""copyright"""

import logging

from utils.utils import common_elements

# iterate through variants found in a proband
# remove any where GQ < 40 if autosome and any DDD_AF > 0.005
# remove any without functional consequence
# return a dict of variants per hgnc_id

def preinheritance_filter(variants):

    variants_per_gene = create_variants_per_gene(variants)
    revel_filter(variants_per_gene)
    return variants_per_gene


def create_variants_per_gene(variants):
    variants_per_gene = {}
    consequences = ['frameshift_variant', 'missense_variant',
                    'splice_donor_variant', 'splice_acceptor_variant',
                    'start_lost', 'stop_gained', 'protein_altering_variant',
                    'transcript_ablation', 'transcript_amplification',
                    'inframe_insertion', 'inframe_deletion', 'stop_lost']

    for v in variants['child'].keys():
        # fail if child GQ < 40
        if int(variants['child'][v].gq) < 40 and variants['child'][
            v].chrom not in ['X', 'Y']:
            logging.info(
                v + " failed low GQ: " + variants['child'][v].gq)
            continue

        # fail if DDD_AF > 0.005
        if float(variants['child'][v].ddd_af) > 0.005:
            logging.info(
                v + " failed high DDD AF: " + variants['child'][v].ddd_af)
            continue

        cqs = variants['child'][v].consequence.split('&')
        coding_cqs = common_elements(cqs, consequences)
        if len(coding_cqs) == 0:
            logging.info(
                v + " failed, no functional consequences: " + variants['child'][
                    v].consequence)
            continue

        # todo fail if the trio genotype is denovo but the variant doesn't pass de novo filters

        hgncid = cqs = variants['child'][v].hgnc_id
        if not hgncid in variants_per_gene.keys():
            variants_per_gene[hgncid] = {}

        variants_per_gene[hgncid][v] = {}
        variants_per_gene[hgncid][v]['child'] = variants['child'][v]
        if v in variants['mum'].keys():
            variants_per_gene[hgncid][v]['mum'] = variants['mum'][v]
        if v in variants['dad'].keys():
            variants_per_gene[hgncid][v]['dad'] = variants['dad'][v]

    return variants_per_gene

def revel_filter(variants_per_gene):
    '''remove missense variants with REVEL < 0.5 unless DNM'''
    for gn in list(variants_per_gene.keys()):
        for varid in list(variants_per_gene[gn].keys()):
            childvar = variants_per_gene[gn][varid]['child']
            if childvar.denovo_snv == "True" or childvar.denovo_indel == "True":
                continue
            elif not childvar.consequence == "missense_variant":
                continue
            elif childvar.revel == '.':
                continue
            else:
                revel = float(childvar.revel)
                if revel < 0.5:
                    logging.info(varid + " failed REVEL filter: " + str(revel))
                    del variants_per_gene[gn][varid]
                    if len(variants_per_gene[gn].keys()) < 1:
                        del variants_per_gene[gn]

