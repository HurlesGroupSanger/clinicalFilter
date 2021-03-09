"""copyright"""

def common_elements(list1, list2):
    return list(set(list1) & set(list2))

def add_single_var_to_candidates(varid, var, hgncid, inh, candidates):
    '''add single variant to candidate variants hash'''
    if not varid in candidates['single_variants'].keys():
        candidates['single_variants'][varid] = {}
        candidates['single_variants'][varid]['mode'] = set()
        candidates['single_variants'][varid]['mode'].add(inh)
        candidates['single_variants'][varid]['variant'] = var
        candidates['single_variants'][varid]['hgncid'] = hgncid
    else:
        candidates['single_variants'][varid]['mode'].add(inh)

def add_compound_het_to_candidates(varid, var, hgncid, inh, candidates):
    '''add compound het to candidate variants hash'''
    # todo could probably streamline this

    if not hgncid in candidates['compound_hets'].keys():
        candidates['compound_hets'][hgncid] = {}
        candidates['compound_hets'][hgncid][varid] = {}
        candidates['compound_hets'][hgncid][varid]['mode'] = set()
        candidates['compound_hets'][hgncid][varid]['mode'].add(inh)
        candidates['compound_hets'][hgncid][varid]['variant'] = var
        candidates['compound_hets'][hgncid][varid]['hgncid'] = hgncid
    else:
        if not varid in candidates['compound_hets'][hgncid].keys():
            candidates['compound_hets'][hgncid][varid] = {}
            candidates['compound_hets'][hgncid][varid]['mode'] = set()
            candidates['compound_hets'][hgncid][varid]['mode'].add(inh)
            candidates['compound_hets'][hgncid][varid]['variant'] = var
            candidates['compound_hets'][hgncid][varid]['hgncid'] = hgncid
        else:
            candidates['compound_hets'][hgncid][varid]['mode'].add(inh)

