import logging

from utils.utils import common_elements
from utils.utils import add_compound_het_to_candidates
from utils.utils import add_single_var_to_candidates

class AllosomalFilter(object):

    def __init__(self, Inheritancefilter, hgncid):
        self.family = Inheritancefilter.family
        self.parents = Inheritancefilter.parents
        self.variants_per_gene = Inheritancefilter.variants_per_gene
        self.candidate_variants = Inheritancefilter.candidate_variants
        self.inheritance_report = Inheritancefilter.inhreport
        self.gene = Inheritancefilter.genes[hgncid]
        self.hgncid = hgncid
        self.proband_X_count = self.family.proband.X_count
        # print(self.family.proband)
        # print(self.family.proband.X_count)
        # exit(0)

    def allosomal_filter(self):
        if self.parents == 'both':
            self.allosomal_both_parents()
        elif self.parents == 'none':
            self.allosomal_no_parents()
        else:
            self.allosomal_single_parent()

    def allosomal_both_parents(self):
        '''screens variants in X/Y where there are both parents and adds candidates
        to candidate_variants
        '''
        variants = self.variants_per_gene[self.hgncid]
        for v in variants.keys():
            # print(v)
            # print(variants[v])
            # print(self.proband_X_count)
            genotype = self.get_variant_genotype(variants[v]['child'], v)
            # if genotype is none then variant fails
            if genotype is None:
                continue
            elif genotype == 'homozygous':
                pass
            elif genotype == 'hemizygous':
                pass
            elif genotype == 'heterozygous':
                pass

    def allosomal_no_parents(self):
        '''screens variants in X/Y where there are no parents and adds candidates
        to candidate_variants
        '''
        variants = self.variants_per_gene[self.hgncid]
        for v in variants.keys():
            genotype = self.get_variant_genotype(variants[v]['child'], v)
            if genotype is None:
                continue
            for inh in self.gene['mode']:
                add_single_var_to_candidates(v, variants[v]['child'], self.hgncid, inh,
                                                 self.candidate_variants)

    def allosomal_single_parent(self):
        '''screens variants in X/Y where there is one parent and adds candidates
        to candidate_variants
        '''
        pass

    def get_variant_genotype(self, variant, v):
        '''work out if a variant is homozygous, heterozygous or hemizygous'''
        if self.proband_X_count >= 2:
            if variant.gt == '0/1':
                genotype = 'heterozygous'
            elif variant.gt == '1/1':
                genotype = 'homozygous'
            else:
                logging.info(v + " fails invalid genotype " + variant.gt)
                return None
        elif self.proband_X_count == 1:
            if variant.gt == '1/1':
                genotype = 'homozygous'
            elif variant.gt == '0/1':
                adsplit = variant.ad.split(',')
                vaf = int(adsplit[1]) / (int(adsplit[0]) + int(adsplit[1]))
                if vaf > 0.9:
                    genotype = 'hemizygous'
                else:
                    logging.info(v + " fails 0/1 variant with low VAF in proband "
                                     "with 1 X chromosome")
                    return None
            else:
                logging.info(v + " fails invalid genotype " + variant.gt)
                return None
        else:
            logging.info(v + " fails invalid X chromsome count")
            return None

        return genotype
