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
            if not variants[v]['child'].is_snv():
                continue
            mum_genotype = variants[v]['child'].get_mum_genotype()
            dad_genotype = variants[v]['child'].get_dad_genotype()
            mum_gt = convert_genotype_to_gt(mum_genotype)
            dad_gt = convert_genotype_to_gt(dad_genotype)
            mum_aff = self.family.mum.affected
            dad_aff = self.family.dad.affected

            genotype = self.get_variant_genotype(variants[v]['child'], v)
            # if dad gt = 0/1 should go to a different (mosaic) pipeline - fail for now
            if dad_gt == '0/1' and variants[v]['child'].chrom == 'X':
                logging.info(v + " failed due to 0/1 paternal genotype in X: ")
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
                    elif inh == 'monoallelic_Y_hem':
                        logging.info(v + "fails inheritance filters homozygous GT in " + inh)
                    else:
                        logging.info(v + " unknown gene mode " + inh)
            elif genotype == 'hemizygous':
                for inh in self.gene['mode']:
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
                    elif inh == 'monoallelic_Y_hem':
                        self.gn_mono_y_hem_gt_hemizygous_parents_filter(
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
                for inh in self.gene['mode']:
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
                    elif inh == 'X-linked over-dominance':
                        self.gn_X_linked_over_dominant_gt_heterozygous_parents_filter(
                            v,
                            variants[
                                v][
                                'child'],
                            mum_gt,
                            dad_gt,
                            mum_aff,
                            dad_aff)
                    elif inh == 'monoallelic_Y_hem':
                        logging.info(v + "fails inheritance filters heterozygous GT in " + inh)
                    else:
                        logging.info(v + " unknown gene mode " + inh)
            else:
                logging.info(v + " unknown genotype " + genotype)

    def allosomal_no_parents(self):
        '''screens variants in X/Y where there are no parents and adds candidates
        to candidate_variants
        '''
        #all allosomal variants on X with no parents pass, Y fail if 0/1
        variants = self.variants_per_gene[self.hgncid]
        for v in variants.keys():
            if not variants[v]['child'].is_snv():
                continue
            genotype = self.get_variant_genotype(variants[v]['child'], v)
            if genotype is None:
                continue
            if variants[v]['child'].chrom == 'X':
                for inh in self.gene['mode']:
                    add_single_var_to_candidates(v, variants[v]['child'], self.hgncid, inh.lower(),
                                                 self.candidate_variants)
            elif variants[v]['child'].chrom == 'Y':
                for inh in self.gene['mode']:
                    if variants[v]['child'].gt == '1/1':
                        add_single_var_to_candidates(v, variants[v]['child'],
                                                 self.hgncid, inh.lower(),
                                                 self.candidate_variants)
                    else:
                        logging.info(
                        v + "fails inheritance filters non-hemizygous GT in " + inh)
            else:
                logging.info(
                    v + "fails inheritance allosomal inheritance filters, chromosome = " + variants[v]['child'].chrom)


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
                                                            'X-linked_dominant',
                                                            'homozygous',
                                                            mum_gt, dad_gt,
                                                            mum_aff, dad_aff)
        logging.info(varid + " failed inheritance filter for homozygous "
                             "variant in X-linked dominant gene")

    def gn_X_linked_over_dominant_gt_homozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        # X-linked over dominant gene, homozygous proband all fail
        self.inheritance_report.populate_inheritance_report('allosomal',
                                                            'X-linked_over_dominant',
                                                            'homozygous',
                                                            mum_gt, dad_gt,
                                                            mum_aff, dad_aff)
        logging.info(varid + " failed inheritance filter for homozygous "
                             "variant in X-linked over dominant gene")

    def gn_hemizygous_gt_hemizygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        #pass unless mum is 1/1 and unaffected
        self.inheritance_report.populate_inheritance_report('allosomal',
                                                            'hemizygous',
                                                            'hemizygous',
                                                            mum_gt, dad_gt,
                                                            mum_aff, dad_aff)
        if not mum_aff and mum_gt == '1/1':
            logging.info(varid + " failed inheritance filter for hemizygous "
                                 "variant in hemizygous gene")
        else:
            add_single_var_to_candidates(varid, var, self.hgncid, 'hemizygous',
                                         self.candidate_variants)

    def gn_X_linked_dominant_gt_hemizygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        # pass unless mum is 1/1 and unaffected
        self.inheritance_report.populate_inheritance_report('allosomal',
                                                            'X-linked_dominant',
                                                            'hemizygous',
                                                            mum_gt, dad_gt,
                                                            mum_aff, dad_aff)
        if not mum_aff and mum_gt == '1/1':
            logging.info(varid + " failed inheritance filter for hemizygous "
                                 "variant in X-linked dominant gene")
        else:
            add_single_var_to_candidates(varid, var, self.hgncid, 'X-linked dominant',
                                         self.candidate_variants)

    def gn_X_linked_over_dominant_gt_hemizygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        # X-linked over dominant gene, hemizygous proband all fail
        self.inheritance_report.populate_inheritance_report('allosomal',
                                                            'X-linked_over_dominant',
                                                            'hemizygous',
                                                            mum_gt, dad_gt,
                                                            mum_aff, dad_aff)
        logging.info(varid + " failed inheritance filter for hemizygous "
                             "variant in X-linked over dominant gene")

    def gn_hemizygous_gt_heterozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        # Hemizygous gene, heterozygous proband. Fail if either parent is 1/1 and unaffected, or if both parents 1/1 and affected
        self.inheritance_report.populate_inheritance_report('allosomal',
                                                            'hemizygous',
                                                            'heterozygous',
                                                            mum_gt, dad_gt,
                                                            mum_aff, dad_aff)
        if not mum_aff and mum_gt == '1/1':
            logging.info(varid + " failed inheritance filter for heterozygous "
                                 "variant in hemizygous gene")
        elif not dad_aff and dad_gt == '1/1':
            logging.info(varid + " failed inheritance filter for heterozygous "
                                 "variant in hemizygous gene")
        elif mum_aff and mum_gt == '1/1' and dad_aff and dad_gt == '1/1':
            logging.info(varid + " failed inheritance filter for heterozygous "
                                 "variant in hemizygous gene")
        else:
            add_single_var_to_candidates(varid, var, self.hgncid, 'hemizygous',
                                         self.candidate_variants)

    def gn_X_linked_dominant_gt_heterozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        # X-linked dominant gene, heterozygous proband. Fail if either parent is 1/1 and unaffected, or if both parents 1/1 and affected
        self.inheritance_report.populate_inheritance_report('allosomal',
                                                            'X-linked_dominant',
                                                            'heterozygous',
                                                             mum_gt, dad_gt,
                                                             mum_aff, dad_aff)
        if not mum_aff and mum_gt == '1/1':
            logging.info(varid + " failed inheritance filter for heterozygous "
                                 "variant in X-linked dominant gene")
        elif not dad_aff and dad_gt == '1/1':
            logging.info(varid + " failed inheritance filter for heterozygous "
                                 "variant in X-linked dominant gene")
        elif mum_aff and mum_gt == '1/1' and dad_aff and dad_gt == '1/1':
            logging.info(varid + " failed inheritance filter for heterozygous "
                                 "variant in X-linked dominant gene")
        else:
            add_single_var_to_candidates(varid, var, self.hgncid, 'X-linked dominant',
                                         self.candidate_variants)


    def gn_X_linked_over_dominant_gt_heterozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        # X-linked over-dominant gene, heterozygous proband. If dad aff pass if DNM
        # if dad unaff pass if DNM or mum_aff and 1/1 or mum unaff and not 0/1
        self.inheritance_report.populate_inheritance_report('allosomal',
                                                            'X-linked_over_dominant',
                                                            'heterozygous',
                                                            mum_gt, dad_gt,
                                                            mum_aff, dad_aff)
        if dad_aff:
            if mum_gt == '0/0' and dad_gt == '0/0':
                add_single_var_to_candidates(varid, var, self.hgncid,
                                             'X-linked over-dominant',
                                             self.candidate_variants)
            else:
                logging.info(varid + " failed inheritance filter for heterozygous "
                                 "variant in X-linked over-dominant gene")
        else:#dad unaff
            if mum_gt == '0/0' and dad_gt == '0/0':
                add_single_var_to_candidates(varid, var, self.hgncid,
                                             'X-linked over-dominant',
                                             self.candidate_variants)
            elif mum_aff and mum_gt == '0/1':
                add_single_var_to_candidates(varid, var, self.hgncid,
                                             'X-linked over-dominant',
                                             self.candidate_variants)
            elif not mum_aff and not mum_gt == '0/1':
                add_single_var_to_candidates(varid, var, self.hgncid,
                                             'X-linked over-dominant',
                                             self.candidate_variants)
            else:
                logging.info(
                    varid + " failed inheritance filter for heterozygous "
                            "variant in X-linked over-dominant gene")

    def gn_mono_y_hem_gt_hemizygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        self.inheritance_report.populate_inheritance_report('allosomal',
                                                            'monoallelic_Y_hemizygous',
                                                            'hemizygous',
                                                            mum_gt, dad_gt,
                                                            mum_aff, dad_aff)
        if not dad_aff and dad_gt == '1/1':
            logging.info(
                varid + " failed inheritance filter for hemizygous "
                        "variant in monoallelic_Y_hemizygous gene")
        else:
            add_single_var_to_candidates(varid, var, self.hgncid,
                                         'monoallelic_Y_hemizygous',
                                         self.candidate_variants)

    def get_variant_genotype(self, variant, v):
        '''work out if a variant is homozygous, heterozygous or hemizygous'''
        if variant.chrom == 'X':
            if self.proband_X_count >= 2:
                if variant.gt == '0/1':
                    genotype = 'heterozygous'
                elif variant.gt == '1/1':
                    genotype = 'hemizygous'
                else:
                    logging.info(v + " fails invalid genotype " + variant.gt)
                    return None
            elif self.proband_X_count == 1:
                if variant.gt == '1/1':
                    genotype = 'hemizygous'
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
        elif variant.chrom == 'Y':
            if variant.gt == '0/1':
                genotype = 'heterozygous'
            elif variant.gt == '1/1':
                genotype = 'hemizygous'
            else:
                logging.info(v + " fails invalid genotype " + variant.gt)
                return None

        return genotype
