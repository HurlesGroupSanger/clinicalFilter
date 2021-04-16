"""copyright"""

import logging

class PostInheritanceFiltering(object):
    '''class for postinheritance filtering'''
    def __init__(self, candidate_variants, family):
        self.candidate_variants = candidate_variants
        self.family = family

    def postinheritance_filter(self):
        '''post-inheritance filtering - MAF and REVEL'''
        self.maf_filter()
        return self.candidate_variants

    def maf_filter(self):
        '''filter non-Biallelic vairants with more stringent MAF thresholds'''
        for v in list(self.candidate_variants['single_variants'].keys()):
            # print(self.candidate_variants['single_variants'][v]['variant'].chrom)
            if self.candidate_variants['single_variants'][v]['mode'] != 'Biallelic':
                ddd_af = self.candidate_variants['single_variants'][v]['variant'].ddd_af
                max_af = self.candidate_variants['single_variants'][v]['variant'].max_af
                if ddd_af == '.':
                    ddd_af = '0'
                if max_af == '.':
                    max_af = '0'
                maximum_af = max(float(ddd_af), float(max_af))

                if self.family.has_both_parents() and maximum_af >= 0.0005:
                    del self.candidate_variants['single_variants'][v]
                    logging.info(
                        v + " failed post-inhertance MAF filter for family with parents, max AF = " + str(
                            maximum_af))
                elif not self.family.has_both_parents() and maximum_af >= 0.0001:
                    del self.candidate_variants['single_variants'][v]
                    logging.info(
                        v + " failed post-inhertance MAF filter for family without parents, max AF = " + str(
                            maximum_af))


