"""copyright"""

import logging


def postinheritance_filter(candidate_variants, family):
    '''post-inheritance filtring - MAF and REVEL'''
    maf_filter(candidate_variants, family)

def maf_filter(candidate_variants, family):
    '''filter non-Biallelic vairants with more stringent MAF thresholds'''
    for v in list(candidate_variants['single_variants'].keys()):
        if candidate_variants['single_variants'][v]['mode'] != 'Biallelic':
            ddd_af = candidate_variants['single_variants'][v]['variant'].ddd_af
            max_af = candidate_variants['single_variants'][v]['variant'].max_af
            if ddd_af == '.':
                ddd_af = '0'
            if max_af == '.':
                max_af = '0'
            maximum_af = max(float(ddd_af), float(max_af))

            if family.has_both_parents() and maximum_af >= 0.0005:
                del candidate_variants['single_variants'][v]
                logging.info(
                    v + " failed post-inhertance MAF filter for family with parents, max AF = " + str(
                        maximum_af))
            elif not family.has_both_parents() and maximum_af >= 0.0001:
                del candidate_variants['single_variants'][v]
                logging.info(
                    v + " failed post-inhertance MAF filter for family without parents, max AF = " + str(
                        maximum_af))


