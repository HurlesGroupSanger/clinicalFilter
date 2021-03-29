import unittest
import tempfile
from family.families import Person, Family
from file_loading.ped_files import openped

class TestLoadPed(unittest.TestCase):

    def setUp(self):
        """ define and write a temporary ped file
        """
        self.maxDiff = None
        self.probandlist = None
        self.tempfile = tempfile.NamedTemporaryFile(mode="w")
        self.path = self.tempfile.name
        self.tempfile.write("fam_ID   proband   dad   mum   XX  2  /path/to/proband_vcf.gz\n")
        self.tempfile.write("fam_ID   dad       0     0     XY  1  /path/to/dad_vcf.gz\n")
        self.tempfile.write("fam_ID   mum       0     0     XX  1  /path/to/mum_vcf.gz\n")
        self.tempfile.flush()

    def test_open_ped_single_family(self):
        """check we can open a ped file with one family
        """
        # load all the components from the file
        self.assertEqual(openped(self.path, self.probandlist), {
            'fam_ID_proband':
                Family(proband=Person(family_id="fam_ID",
                                      person_id="proband",
                                      dad_id="dad", mum_id="mum",
                                      sex="XX", affected="2",
                                      path="/path/to/proband_vcf.gz"),
                       dad=Person(family_id="fam_ID",
                                  person_id="dad", dad_id="0",
                                  mum_id="0", sex="XY", affected="1",
                                  path="/path/to/dad_vcf.gz"),
                       mum=Person(family_id="fam_ID",
                                  person_id="mum", dad_id="0",
                                  mum_id="0", sex="XX", affected="1",
                                  path="/path/to/mum_vcf.gz"))})

    def test_open_ped_multiple_sibs(self):
        """ check that we correctly parse a ped file with multiple siblings
        """
        # add an extra sibling
        self.tempfile.write("fam_ID  sib   dad  mum  XX  2  /path/to/sib_vcf.gz\n")
        self.tempfile.flush()

        # load all the components from the file
        self.assertEqual(openped(self.path, self.probandlist), {
            'fam_ID_proband':
                Family(proband=Person(family_id="fam_ID",
                                      person_id="proband",
                                      dad_id="dad", mum_id="mum",
                                      sex="XX", affected="2",
                                      path="/path/to/proband_vcf.gz"),
                       dad=Person(family_id="fam_ID",
                                  person_id="dad", dad_id="0",
                                  mum_id="0", sex="XY", affected="1",
                                  path="/path/to/dad_vcf.gz"),
                       mum=Person(family_id="fam_ID",
                                  person_id="mum", dad_id="0",
                                  mum_id="0", sex="XX", affected="1",
                                  path="/path/to/mum_vcf.gz")),
            'fam_ID_sib': Family(proband=Person(family_id="fam_ID",
                                                person_id="sib",
                                                dad_id="dad",
                                                mum_id="mum",
                                                sex="XX",
                                                affected="2",
                                                path="/path/to/sib_vcf.gz"),
                                 dad=Person(family_id="fam_ID",
                                            person_id="dad",
                                            dad_id="0",
                                            mum_id="0", sex="XY",
                                            affected="1",
                                            path="/path/to/dad_vcf.gz"),
                                 mum=Person(family_id="fam_ID",
                                            person_id="mum",
                                            dad_id="0",
                                            mum_id="0", sex="XX",
                                            affected="1",
                                            path="/path/to/mum_vcf.gz"))})


    def test_open_ped_multiple_families(self):
        """ check that we correctly parse a ped file with multiple families
        """
        # add an extra family, with multiple sibs
        self.tempfile.write("fam_ID2  proband2  dad2  mum2  XX  2  /path/to/proband2_vcf.gz\n")
        self.tempfile.write("fam_ID2  dad2      0     0     XY  1  /path/to/dad2_vcf.gz\n")
        self.tempfile.write("fam_ID2  mum2      0     0     XX  1  /path/to/mum2_vcf.gz\n")
        self.tempfile.write("fam_ID2  sib       dad2  mum2  XX  2  /path/to/sib_vcf.gz\n")
        self.tempfile.flush()

        # load all the components from the file
        self.assertEqual(openped(self.path, self.probandlist), {
            'fam_ID_proband':
                Family(proband=Person(family_id="fam_ID",
                                      person_id="proband",
                                      dad_id="dad", mum_id="mum",
                                      sex="XX", affected="2",
                                      path="/path/to/proband_vcf.gz"),
                       dad=Person(family_id="fam_ID",
                                  person_id="dad", dad_id="0",
                                  mum_id="0", sex="XY", affected="1",
                                  path="/path/to/dad_vcf.gz"),
                       mum=Person(family_id="fam_ID",
                                  person_id="mum", dad_id="0",
                                  mum_id="0", sex="XX", affected="1",
                                  path="/path/to/mum_vcf.gz")),
            'fam_ID2_proband2': Family(proband=Person(family_id="fam_ID2",
                                                      person_id="proband2",
                                                      dad_id="dad2", mum_id="mum2",
                                                      sex="XX", affected="2",
                                                      path="/path/to/proband2_vcf.gz"),
                                       dad=Person(family_id="fam_ID2",
                                                  person_id="dad2", dad_id="0",
                                                  mum_id="0", sex="XY", affected="1",
                                                  path="/path/to/dad2_vcf.gz"),
                                       mum=Person(family_id="fam_ID2",
                                                  person_id="mum2", dad_id="0",
                                                  mum_id="0", sex="XX", affected="1",
                                                  path="/path/to/mum2_vcf.gz")),
            'fam_ID2_sib': Family(proband=Person(family_id="fam_ID2",
                                                 person_id="sib",
                                                 dad_id="dad2", mum_id="mum2",
                                                 sex="XX", affected="2",
                                                 path="/path/to/sib_vcf.gz"),
                                  dad=Person(family_id="fam_ID2",
                                             person_id="dad2", dad_id="0",
                                             mum_id="0", sex="XY", affected="1",
                                             path="/path/to/dad2_vcf.gz"),
                                  mum=Person(family_id="fam_ID2",
                                             person_id="mum2", dad_id="0",
                                             mum_id="0", sex="XX", affected="1",
                                             path="/path/to/mum2_vcf.gz"))
        })

    def test_families_from_list(self):
        """check families created from a proband list"""
        # add an extra family, with multiple sibs
        self.tempfile.write("fam_ID2  proband2  dad2  mum2  XX  2  /path/to/proband2_vcf.gz\n")
        self.tempfile.write("fam_ID2  dad2      0     0     XY  1  /path/to/dad2_vcf.gz\n")
        self.tempfile.write("fam_ID2  mum2      0     0     XX  1  /path/to/mum2_vcf.gz\n")
        self.tempfile.write("fam_ID2  sib       dad2  mum2  XX  2  /path/to/sib_vcf.gz\n")
        self.tempfile.flush()
        #add the list of probands
        self.tempfile2 = tempfile.NamedTemporaryFile(mode="w")
        self.path2 = self.tempfile2.name
        self.tempfile2.write("proband\n")
        self.tempfile2.write("proband2\n")
        self.tempfile2.write("sib\n")
        self.tempfile2.flush()

        self.assertEqual(openped(self.path, self.probandlist), {
            'fam_ID_proband':
                Family(proband=Person(family_id="fam_ID",
                                      person_id="proband",
                                      dad_id="dad", mum_id="mum",
                                      sex="XX", affected="2",
                                      path="/path/to/proband_vcf.gz"),
                       dad=Person(family_id="fam_ID",
                                  person_id="dad", dad_id="0",
                                  mum_id="0", sex="XY", affected="1",
                                  path="/path/to/dad_vcf.gz"),
                       mum=Person(family_id="fam_ID",
                                  person_id="mum", dad_id="0",
                                  mum_id="0", sex="XX", affected="1",
                                  path="/path/to/mum_vcf.gz")),
            'fam_ID2_proband2': Family(proband=Person(family_id="fam_ID2",
                                                      person_id="proband2",
                                                      dad_id="dad2", mum_id="mum2",
                                                      sex="XX", affected="2",
                                                      path="/path/to/proband2_vcf.gz"),
                                       dad=Person(family_id="fam_ID2",
                                                  person_id="dad2", dad_id="0",
                                                  mum_id="0", sex="XY", affected="1",
                                                  path="/path/to/dad2_vcf.gz"),
                                       mum=Person(family_id="fam_ID2",
                                                  person_id="mum2", dad_id="0",
                                                  mum_id="0", sex="XX", affected="1",
                                                  path="/path/to/mum2_vcf.gz")),
            'fam_ID2_sib': Family(proband=Person(family_id="fam_ID2",
                                                 person_id="sib",
                                                 dad_id="dad2", mum_id="mum2",
                                                 sex="XX", affected="2",
                                                 path="/path/to/sib_vcf.gz"),
                                  dad=Person(family_id="fam_ID2",
                                             person_id="dad2", dad_id="0",
                                             mum_id="0", sex="XY", affected="1",
                                             path="/path/to/dad2_vcf.gz"),
                                  mum=Person(family_id="fam_ID2",
                                             person_id="mum2", dad_id="0",
                                             mum_id="0", sex="XX", affected="1",
                                             path="/path/to/mum2_vcf.gz"))
        })

    def test_family_parents(self):
        """check the functions to determine what parents a family has"""
        testproband1 = Person("fam_ID", "proband", "0", "0", "XX", "2",
                           "/path/to/proband_vcf.gz")
        testproband2 = Person("fam_ID", "proband", "dad", "mum", "XX", "2",
                           "/path/to/proband_vcf.gz")
        testproband3 = Person("fam_ID", "proband", "0", "mum", "XX", "2",
                           "/path/to/proband_vcf.gz")
        testproband4 = Person("fam_ID", "proband", "dad", "0", "XX", "2",
                              "/path/to/proband_vcf.gz")
        testmum = Person("fam_ID", "mum", "0", "0", "XX", "2",
                              "/path/to/mum_vcf.gz")
        testdad = Person("fam_ID", "dad", "0", "0", "XY", "2",
                         "/path/to/dad_vcf.gz")
        testfamily1 = Family(testproband1, None, None)
        testfamily2 = Family(testproband2, testmum, testdad)
        testfamily3 = Family(testproband3, testmum, None)
        testfamily4 = Family(testproband4, None, testdad)
        self.assertTrue(testfamily1.has_no_parents())
        self.assertFalse(testfamily1.has_both_parents())
        self.assertFalse(testfamily1.has_mum())
        self.assertFalse(testfamily1.has_dad())
        self.assertTrue(testfamily2.has_both_parents())
        self.assertTrue(testfamily2.has_mum())
        self.assertTrue(testfamily2.has_dad())
        self.assertFalse(testfamily2.has_no_parents())
        self.assertFalse(testfamily3.has_both_parents())
        self.assertTrue(testfamily3.has_mum())
        self.assertFalse(testfamily3.has_dad())
        self.assertFalse(testfamily3.has_no_parents())
        self.assertFalse(testfamily4.has_both_parents())
        self.assertFalse(testfamily4.has_mum())
        self.assertTrue(testfamily4.has_dad())
        self.assertFalse(testfamily4.has_no_parents())


    def test_affected_status(self):
        """check that affected status is parsed correctly in person object
        """
        affperson = Person("fam_ID", "proband", "dad", "mum", "XX", "2",
                           "/path/to/proband_vcf.gz")
        unaffperson = Person("fam_ID", "unaff", "dad", "mum", "XX", "1",
                           "/path/to/unaff_vcf.gz")
        affstatus = affperson.get_affected_status()
        unaffstatus = unaffperson.get_affected_status()
        self.assertTrue(affstatus)
        self.assertFalse(unaffstatus)

if __name__ == '__main__':
    unittest.main()
