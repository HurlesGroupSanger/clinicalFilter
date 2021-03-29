"""copyright"""

import logging
from family.families import Family, Person

def openped(filename, proband_list):
    """
    read ped files and return a dict of family objects
    """
    families = {}
    people = {}
    parents = {}
    probands = {}
    with open(filename, 'r') as p:
        lines = p.readlines()
        for l in lines:
            #create a person object for each individual
            individual = Person(*l.strip().split())
            #populate dict of person id:object
            person_id = individual.get_id()
            people[person_id] = individual
            #add parents to the dict of parents in this ped file
            parent_ids = individual.get_parents()
            for parent_id in parent_ids:
                parents[parent_id] = 1

    #get a list of probands to analyse and create family objects
    if proband_list:
        """parse probands from a list"""
        with open(proband_list, 'r') as p:
            lines = p.readlines()
            for l in lines:
                person = l.strip()
                if person in people:
                    probands[person] = 1
                else:
                    logging.debug(person + " not found in ped file, skipping")
    else:
        """assume we have only single generation families - so the probands are
         the affected individuals who are not parents
         """
        for person in people.keys():
            affected = people[person].get_affected_status()
            if affected and person not in parents.keys():
                probands[person] = 1

    for probandid in probands.keys():
        """create families. Family id plus proband id should be unique so 
        combine this for a key"""
        if probandid in people.keys():
            famid = people[probandid].get_family_id() + "_" + probandid
            mumid = people[probandid].get_mum_id()
            dadid = people[probandid].get_dad_id()

            if mumid == '0':
                mum = None
            elif mumid in people.keys():
                mum = people[mumid]
            else:
                mum = None
                logging.debug("Proband " + probandid + ": mother not found in ped "
                                "file - skipping mother " + mumid)
            if dadid == '0':
                dad = None
            elif dadid in people.keys():
                dad = people[dadid]
            else:
                dad = None
                logging.debug("Proband " + probandid + ": father not found in ped "
                                "file - skipping father " + dadid)

            family = Family(people[probandid], mum, dad)
            families[famid] = family
        else:
            logging.debug("Proband " + probandid + " not found in ped file - skipping")

    return families

def create_ped(pedfile, child, mother, father, sex, mum_aff, dad_aff):
    """create ped file where there isn't one inputted"""
    with open(pedfile, 'w') as p:
        mum = ''
        dad = ''
        if mother is None:
            mum = "0"
        else:
            mum = "mum"
        if father is None:
            dad = "0"
        else:
            dad = "dad"
        probandline = ("\t").join(["family", "child", dad, mum, sex, "2", child]) + "\n"
        p.write(probandline)
        if father:
            dadline = ("\t").join(["family", "dad", "0", "0", "XY", dad_aff, father]) + "\n"
            p.write(dadline)
        if mother:
            mumline = ("\t").join(["family", "mum", "0", "0", "XX", mum_aff, mother]) + "\n"
            p.write(mumline)