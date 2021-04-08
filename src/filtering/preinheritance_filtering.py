"""copyright"""

import logging

from utils.utils import common_elements

class PreInheritanceFiltering(object):
    '''class for preinheritance filtering'''

    # iterate through variants found in a proband
    # remove any where GQ < 40 if autosome and any DDD_AF > 0.005
    # remove any without functional consequence
    # return a dict of variants per hgnc_id

    def __init__(self, variants):
        self.variants = variants

    def preinheritance_filter(self):

        variants_per_gene = self.create_variants_per_gene()
        self.revel_filter(variants_per_gene)
        return variants_per_gene


    def create_variants_per_gene(self):
        variants_per_gene = {}
        consequences = ['frameshift_variant', 'missense_variant',
                        'splice_donor_variant', 'splice_acceptor_variant',
                        'start_lost', 'stop_gained', 'protein_altering_variant',
                        'transcript_ablation', 'transcript_amplification',
                        'inframe_insertion', 'inframe_deletion', 'stop_lost']

        for v in self.variants['child'].keys():
            # fail if child GQ < 40
            if int(self.variants['child'][v].gq) < 40 and self.variants['child'][
                v].chrom not in ['X', 'Y']:
                logging.info(
                    v + " failed low GQ: " + self.variants['child'][v].gq)
                continue

            # fail if DDD_AF > 0.005
            if float(self.variants['child'][v].ddd_af) > 0.005:
                logging.info(
                    v + " failed high DDD AF: " + self.variants['child'][v].ddd_af)
                continue

            cqs = self.variants['child'][v].consequence.split('&')
            coding_cqs = common_elements(cqs, consequences)
            if len(coding_cqs) == 0:
                logging.info(
                    v + " failed, no functional consequences: " + self.variants['child'][
                        v].consequence)
                continue

            # todo fail if the trio genotype is denovo but the variant doesn't pass de novo filters

            hgncid = cqs = self.variants['child'][v].hgnc_id
            if not hgncid in variants_per_gene.keys():
                variants_per_gene[hgncid] = {}

            variants_per_gene[hgncid][v] = {}
            variants_per_gene[hgncid][v]['child'] = self.variants['child'][v]
            if v in self.variants['mum'].keys():
                variants_per_gene[hgncid][v]['mum'] = self.variants['mum'][v]
            if v in self.variants['dad'].keys():
                variants_per_gene[hgncid][v]['dad'] = self.variants['dad'][v]

        return variants_per_gene

    def revel_filter(self, variants_per_gene):
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

