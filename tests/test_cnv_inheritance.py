import unittest
import copy

from tests.test_utils import create_test_candidate_vars
from tests.test_utils import create_test_person
from tests.test_utils import create_test_family
from tests.test_utils import create_test_variants_per_gene

class TestAutosomalInheritanceFilter(unittest.TestCase):

    def setUp(self):
        self.maxDiff = None

    def test_non_ddg2p_filter(self):
        #should pass if cnv length > 1000000, fail if smaller and no ddg2p overlap
        pass

