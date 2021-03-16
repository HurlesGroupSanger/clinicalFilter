"""copyright"""

# methods to create person, family and variant objects to use in tests

from variants.snv import SNV
from family.families import Person
from family.families import Family

def create_test_person(family_id, person_id, dad_id, mum_id, sex, affected, path):
    person = Person(family_id, person_id, dad_id, mum_id, sex, affected, path)
    return person

def create_test_family(child, mum, dad):
    family = Family(child, mum, dad)
    return family

def create_test_snv(vardata):
    '''create an SNV from a variant hash'''
    Var = SNV
    var = Var(vardata)
    return var

def create_test_candidate_vars(single_vars, compound_hets):
    '''create candidate variants hash from variant data'''
    candidates = {'single_variants': {}, 'compound_hets': {}}

    for vid in single_vars.keys():
        Var = SNV
        var = Var(single_vars[vid]['variant'])
        candidates['single_variants'][vid] = {}
        candidates['single_variants'][vid]['variant'] = var
        candidates['single_variants'][vid]['mode'] = single_vars[vid]['mode']

    for gn in compound_hets.keys():
        candidates['compound_hets'][gn] = {}
        for cvid in compound_hets[gn].keys():
            Var = SNV
            var = Var(compound_hets[gn][cvid]['variant'])
            candidates['compound_hets'][gn][cvid]['variant'] = var
            candidates['compound_hets'][gn][cvid]['mode'] = compound_hets[gn][cvid]['mode']

    return candidates
