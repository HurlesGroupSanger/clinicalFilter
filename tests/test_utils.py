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