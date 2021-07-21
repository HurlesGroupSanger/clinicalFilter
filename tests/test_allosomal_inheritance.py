import unittest

from tests.test_utils import create_test_candidate_vars
from tests.test_utils import create_test_person
from tests.test_utils import create_test_family
from tests.test_utils import create_test_variants_per_gene

from filtering.inheritance_filtering import InheritanceFiltering


class TestAllosomalInheritanceFilter(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None

    def test_allosomal_no_parents(self):
        pass

    def test_get_variant_genotype(self):
        #test that variant genotype is calculated correctly from child chromosomal sex
        pass

    def test_gn_hemizygous_gt_homozygous_parents_filter(self):
        pass

    def test_gn_X_linked_dominant_gt_homozygous_parents_filter(self):
        pass

    def test_gn_X_linked_over_dominant_gt_homozygous_parents_filter(self):
        pass

    def test_gn_hemizygous_gt_hemizygous_parents_filter(self):
        pass

    def test_gn_X_linked_dominant_gt_hemizygous_parents_filter(self):
        pass

    def test_gn_X_linked_over_dominant_gt_hemizygous_parents_filter(self):
        pass

    def test_gn_hemizygous_gt_heterozygous_parents_filter(self):
        pass

    def test_gn_X_linked_dominant_gt_heterozygous_parents_filter(self):
        pass

    def test_gn_X_linked_over_dominant_gt_heterozygous_parents_filter(self):
        pass