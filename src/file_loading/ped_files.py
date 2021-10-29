"""
Copyright (c) 2021 Genome Research Limited
Author: Ruth Eberhardt <re3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

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
            # create a person object for each individual
            individual = Person(*l.strip().split())
            # populate dict of person id:object
            person_id = individual.get_id()
            people[person_id] = individual
            # add parents to the dict of parents in this ped file
            parent_ids = individual.get_parents()
            for parent_id in parent_ids:
                parents[parent_id] = 1

    # get a list of probands to analyse and create family objects
    if proband_list:
        # parse probands from a list
        with open(proband_list, 'r') as p:
            lines = p.readlines()
            for l in lines:
                person = l.strip()
                if person in people:
                    probands[person] = 1
                else:
                    logging.debug(person + " not found in ped file, skipping")
    else:
        # assume we have only single generation families - so the probands are
        # the affected individuals who are not parents
        for person in people.keys():
            affected = people[person].get_affected_status()
            if affected and person not in parents.keys():
                probands[person] = 1

    for probandid in probands.keys():
        # create families. Family id plus proband id should be unique so
        # combine this for a key
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
                logging.debug("Proband " + probandid
                              + ": mother not found in ped "
                                "file - skipping mother " + mumid)
            if dadid == '0':
                dad = None
            elif dadid in people.keys():
                dad = people[dadid]
            else:
                dad = None
                logging.debug("Proband " + probandid
                              + ": father not found in ped "
                                "file - skipping father " + dadid)

            family = Family(people[probandid], mum, dad)
            families[famid] = family
        else:
            logging.debug("Proband " + probandid
                          + " not found in ped file - skipping")

    return families

def create_ped(pedfile, child, mother, father, sex, mum_aff, dad_aff):
    """
    create ped file where there isn't one inputted
    """
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