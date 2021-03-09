"""
copyright
"""
import logging

class Person(object):
    """
    Person object: ID, VCF, gender, affected status
    """
    def __init__(self, family_id, person_id, dad_id, mum_id, sex, affected, path):
        self.family_id = family_id#note that a person may be in >1 family if it is a parent
        self.person_id = person_id
        self.mum_id = mum_id
        self.dad_id = dad_id
        self.vcf_path = path
        self.gender = sex
        #convert affected to true/false
        if affected == '2':
            self.affected = True
        elif affected == '1':
            self.affected = False
        else:
            raise ValueError("Unknown affected status: " + affected + " should be '1' or '2'")

    def __repr__(self):
        return 'Person(person_id="{}", dad_id="{}", mum_id="{}", gender="{}", ' \
               'affected="{}", path="{}")'.format(self.get_id(), \
                self.get_dad_id(), self.get_mum_id(), self.get_gender(), \
                self.get_affected_status(), self.get_vcf_path())

    def get_id(self):
        """
        get person id
        """
        return self.person_id

    def get_mum_id(self):
        """
        get person id
        """
        return self.mum_id

    def get_dad_id(self):
        """
        get person id
        """
        return self.dad_id

    def get_family_id(self):
        """
        get family id
        """
        return self.family_id

    def get_gender(self):
        """
        get gender
        """
        return self.gender

    def get_affected_status(self):
        """
        get affected status
        """
        return self.affected

    def get_vcf_path(self):
        """
        get path to vcf file for individual
        """
        return self.vcf_path

    def get_parents(self):
        """get ids of parents (if any)
        """
        parent_ids = []
        if not self.mum_id == "0":
            parent_ids.append(self.mum_id)
        if not self.dad_id == "0":
            parent_ids.append(self.dad_id)
        return parent_ids

    def __eq__(self, other):
        return self.family_id == other.family_id and \
               self.person_id == other.person_id and \
               self.mum_id == other.mum_id and \
               self.dad_id == other.dad_id and \
               self.affected == other.affected and \
               self.vcf_path == other.vcf_path and \
               self.gender == other.gender



class Family(object):
    """
    Family object: Child, Mum, Dad person objects. Can also be proband only or
    single parent
    """

    def __init__(self, proband, mum, dad):
        self.proband = proband
        self.mum = mum
        self.dad = dad

    def __repr__(self):
        return 'Family(proband={}, dad={}, mum={})'.format(self.proband, self.dad, \
                self.mum)

    def has_mum(self):
        """do we have the mum? returns true/false"""
        if self.mum is None:
            return False
        else:
            return True

    def has_dad(self):
        """do we have the dad? returns true/false"""
        if self.dad is None:
            return False
        else:
            return True

    def has_both_parents(self):
        """do we have both parents? returns true/false"""
        if self.mum is None or self.dad is None:
            return False
        else:
            return True

    def has_no_parents(self):
        """do we have neither parent? returns true/false"""
        if self.dad is None and self.mum is None:
            return True
        else:
            return False

    def __eq__(self, other):
        return self.proband == other.proband and \
               self.mum == other.mum and \
               self.dad == other.dad


