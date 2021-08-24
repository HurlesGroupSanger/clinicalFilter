"""copyright"""

import logging
from itertools import combinations


class CompoundHetScreen(object):

    def __init__(self, candidate_variants, family):
        self.candidate_variants = candidate_variants
        self.family = family

    def screen_compound_hets(self):
        # sort out compound hets
        compound_het_passes = {}
        for gn in self.candidate_variants['compound_hets'].keys():
            # print(gn)
            if len(self.candidate_variants['compound_hets'][gn].keys()) < 2:
                for v in self.candidate_variants['compound_hets'][gn].keys():
                    logging.info(
                        v + " failed compound het screen: <2 vars in hgnc " + gn)

            combs_to_screen = list(
                combinations(
                    self.candidate_variants['compound_hets'][gn].keys(), 2))
            for pair in combs_to_screen:
                var1 = self.candidate_variants['compound_hets'][gn][pair[0]][
                    'variant']
                var2 = self.candidate_variants['compound_hets'][gn][pair[1]][
                    'variant']

                if var1 != var2:
                    if self.is_compound_pair(pair[0], var1, pair[1], var2):
                        if not gn in compound_het_passes.keys():
                            compound_het_passes[gn] = {}
                        if not pair[0] in compound_het_passes[gn].keys():
                            compound_het_passes[gn][pair[0]] = {}
                        if not pair[1] in compound_het_passes[gn].keys():
                            compound_het_passes[gn][pair[1]] = {}

                        compound_het_passes[gn][pair[0]]['variant'] = var1
                        compound_het_passes[gn][pair[0]]['mode'] = \
                            self.candidate_variants['compound_hets'][gn][
                                pair[0]][
                                'mode']
                        compound_het_passes[gn][pair[1]]['variant'] = var2
                        compound_het_passes[gn][pair[1]]['mode'] = \
                            self.candidate_variants['compound_hets'][gn][
                                pair[1]][
                                'mode']

        self.candidate_variants['compound_hets'] = compound_het_passes
        # for gn in candidate_variants['compound_hets'].keys():
        #     print(gn)
        # print(candidate_variants['compound_hets'])
        # return self.candidate_variants

    def is_compound_pair(self, varid1, var1, varid2, var2):
        '''Test to see if a pair of variants could be a compound het'''

        if self.family.has_no_parents():
            """If there are no parents and the variants are not both missense then 
            pass"""
            if var1.consequence == 'missense_variant' and \
                    var2.consequence == 'missense_variant':
                logging.info(varid1 + " " + varid2 + " failed compound het "
                                                     "screen, no parents and both missense")
                return False
            elif (var1.pid != '.' and var2.pid != '.') and (
                    var1.pid == var2.pid):
                logging.info(varid1 + " " + varid2 + " failed compound het "
                                                     "screen, variants in cis")
                return False

            else:
                return True
        elif self.family.has_both_parents():
            if var1.chrom == 'X' and not self.family.dad.get_affected_status() and \
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
            elif (var1.dnm == "DNM" and \
                  var2.get_mum_genotype() != '0' and \
                  var2.get_dad_genotype() == '0') or \
                    (var1.dnm == "DNM" and \
                     var2.get_mum_genotype() == '0' and \
                     var2.get_dad_genotype() != '0') or \
                    (var2.dnm == "DNM" and \
                     var1.get_mum_genotype() != '0' and \
                     var1.get_dad_genotype() == '0') or \
                    (var2.dnm == "DNM" and \
                     var1.get_mum_genotype() == '0' and \
                     var1.get_dad_genotype() != '0'):
                '''one variant is DNM and the other is inherited'''
                return True
            else:
                logging.info(
                    varid1 + " " + varid2 + " failed compound het screen")
                return False

        elif self.family.has_dad():
            pass
        # todo screening of compound hets with one parent
        elif self.family.has_mum():
            pass
        else:
            print(
                "Error: should not get here as family must either have both, one \
                 or no parents")
            logging.error("Can't parse compound hets - family error")
            exit(1)
        logging.info(varid1 + " " + varid2 + " failed compound het screen")
        return False
