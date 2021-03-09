"""copyright"""

import logging


def postinheritance_filter(candidate_variants, family):
    '''post-inheritance filtring - MAF and REVEL'''
    maf_filter(candidate_variants, family)
    # revel_filter(candidate_variants)
    # #further check for compound hets
    # final_compound_het_filter(candidate_variants, family)


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


def revel_filter(candidate_variants):
    '''filter on REVEL prediction'''
    for v in list(candidate_variants['single_variants'].keys()):
        if candidate_variants['single_variants'][v][
            'variant'].denovo_snv == "True" or \
                candidate_variants['single_variants'][v][
                    'variant'].denovo_indel == "True":
            # DNMs are not filtered on REVEL
            continue
        elif not candidate_variants['single_variants'][v]['variant'].consequence == 'missense_variant':
            # only missense variants are filtered with polyphen
            continue
        elif candidate_variants['single_variants'][v]['variant'].revel == '.':
            # ignore those with no revel annotation
            continue
        else:
            revel = float(candidate_variants['single_variants'][v]['variant'].revel)
            if revel < 0.5:
                del candidate_variants['single_variants'][v]
                logging.info(v + " failed post-inhertance flter, REVEL = " + str(
                        revel))
    for gn in list(candidate_variants['compound_hets'].keys()):
        for v in list(candidate_variants['compound_hets'][gn].keys()):
            if candidate_variants['compound_hets'][gn][v]['variant'].denovo_snv == "True" or candidate_variants['compound_hets'][gn][v]['variant'].denovo_indel == "True":
                # DNMs are not filtered on REVEL
                continue
            elif not candidate_variants['compound_hets'][gn][v]['variant'].consequence == 'missense_variant':
                # only missense variants are filtered with polyphen
                continue
            elif candidate_variants['compound_hets'][gn][v]['variant'].revel == '.':
                # ignore those with no revel annotation
                continue
            else:
                revel = float(candidate_variants['compound_hets'][gn][v]['variant'].revel)
                if revel < 0.5:
                    del candidate_variants['compound_hets'][gn][v]
                    logging.info(
                        v + " failed post-inhertance flter, REVEL = " + str(
                            revel))
        # if there are <2 variants left in a pair then reject
        numvars_gene = len(candidate_variants['compound_hets'][gn].keys())
        if numvars_gene < 2:
            for vid in candidate_variants['compound_hets'][gn].keys():
                logging.info(vid + " failed post-inhertance flter, compound het partner failed REVEL")
            del candidate_variants['compound_hets'][gn]


# def final_compound_het_filter(candidate_variants, family):
#     '''final compound het filter, now variants have been removed by REVEL
#     filter are all compound hets still proper pairs'''
#     pass
