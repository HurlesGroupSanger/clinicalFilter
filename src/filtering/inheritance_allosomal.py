import logging

from utils.utils import add_single_var_to_candidates
from utils.utils import convert_genotype_to_gt

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
            mum_genotype = variants[v]['child'].get_mum_genotype()
            dad_genotype = variants[v]['child'].get_dad_genotype()
            mum_gt = convert_genotype_to_gt(mum_genotype)
            dad_gt = convert_genotype_to_gt(dad_genotype)
            mum_aff = self.family.mum.affected
            dad_aff = self.family.dad.affected

            genotype = self.get_variant_genotype(variants[v]['child'], v)
            # if genotype is none then variant fails
            if genotype is None:
                continue
            elif genotype == 'homozygous':
                for inh in self.gene['mode']:
                    if inh == 'Hemizygous':
                        self.gn_hemizygous_gt_homozygous_parents_filter(v, variants[v]['child'], mum_gt, dad_gt, mum_aff, dad_aff)
                    elif inh == 'X-linked dominant':
                        self.gn_X_linked_dominant_gt_homozygous_parents_filter(v,
                                                                        variants[
                                                                            v][
                                                                            'child'],
                                                                        mum_gt,
                                                                        dad_gt,
                                                                        mum_aff,
                                                                        dad_aff)
                    elif inh == 'X-linked over-dominant':
                        self.gn_X_linked_over_dominant_gt_homozygous_parents_filter(v,
                                                                        variants[
                                                                            v][
                                                                            'child'],
                                                                        mum_gt,
                                                                        dad_gt,
                                                                        mum_aff,
                                                                        dad_aff)
                    else:
                        logging.info(v + " unknown gene mode " + inh)
            elif genotype == 'hemizygous':
                if inh == 'Hemizygous':
                    self.gn_hemizygous_gt_hemizygous_parents_filter(v,
                                                                    variants[v][
                                                                        'child'],
                                                                    mum_gt,
                                                                    dad_gt,
                                                                    mum_aff,
                                                                    dad_aff)
                elif inh == 'X-linked dominant':
                    self.gn_X_linked_dominant_gt_hemizygous_parents_filter(v,
                                                                           variants[
                                                                               v][
                                                                               'child'],
                                                                           mum_gt,
                                                                           dad_gt,
                                                                           mum_aff,
                                                                           dad_aff)
                elif inh == 'X-linked over-dominant':
                    self.gn_X_linked_over_dominant_gt_hemizygous_parents_filter(
                        v,
                        variants[
                            v][
                            'child'],
                        mum_gt,
                        dad_gt,
                        mum_aff,
                        dad_aff)
                else:
                    logging.info(v + " unknown gene mode " + inh)
            elif genotype == 'heterozygous':
                if inh == 'Hemizygous':
                    self.gn_hemizygous_gt_heterozygous_parents_filter(v,
                                                                    variants[v][
                                                                        'child'],
                                                                    mum_gt,
                                                                    dad_gt,
                                                                    mum_aff,
                                                                    dad_aff)
                elif inh == 'X-linked dominant':
                    self.gn_X_linked_dominant_gt_heterozygous_parents_filter(v,
                                                                           variants[
                                                                               v][
                                                                               'child'],
                                                                           mum_gt,
                                                                           dad_gt,
                                                                           mum_aff,
                                                                           dad_aff)
                elif inh == 'X-linked over-dominant':
                    self.gn_X_linked_over_dominant_gt_heterozygous_parents_filter(
                        v,
                        variants[
                            v][
                            'child'],
                        mum_gt,
                        dad_gt,
                        mum_aff,
                        dad_aff)
                else:
                    logging.info(v + " unknown gene mode " + inh)
            else:
                logging.info(v + " unknown genotype " + genotype)

    def allosomal_no_parents(self):
        '''screens variants in X/Y where there are no parents and adds candidates
        to candidate_variants
        '''
        #all allosomal variants with no parents pass
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

    def gn_hemizygous_gt_homozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        #hemizygous gene, homozygous proband all fail
        self.inheritance_report.populate_inheritance_report('allosomal',
                                                            'hemizygous', 'homozygous',
                                                            mum_gt, dad_gt,
                                                            mum_aff, dad_aff)
        logging.info(varid + " failed inheritance filter for homozygous "
                             "variant in hemizygous gene")

    def gn_X_linked_dominant_gt_homozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        # X-linked dominant gene, homozygous proband all fail
        self.inheritance_report.populate_inheritance_report('allosomal',
                                                            'X-linked dominant',
                                                            'homozygous',
                                                            mum_gt, dad_gt,
                                                            mum_aff, dad_aff)
        logging.info(varid + " failed inheritance filter for homozygous "
                             "variant in X-linked dominant gene")

    def gn_X_linked_over_dominant_gt_homozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        # X-linked over dominant gene, homozygous proband all fail
        self.inheritance_report.populate_inheritance_report('allosomal',
                                                            'X-linked over-dominant',
                                                            'homozygous',
                                                            mum_gt, dad_gt,
                                                            mum_aff, dad_aff)
        logging.info(varid + " failed inheritance filter for homozygous "
                             "variant in X-linked over dominant gene")

    def gn_hemizygous_gt_hemizygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        pass

    def gn_X_linked_dominant_gt_hemizygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        pass

    def gn_X_linked_over_dominant_gt_hemizygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        pass

    def gn_hemizygous_gt_heterozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        pass

    def gn_X_linked_dominant_gt_heterozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        pass

    def gn_X_linked_over_dominant_gt_heterozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
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
