import unittest
import tempfile

from file_loading.load_genes_and_regions import load_genes


class TestLoadGenesRegions(unittest.TestCase):

    def setUp(self):
        '''create genes file to load'''
        self.maxDiff = None
        header = ("\t").join(
            ['chr', 'start', 'stop', 'gene', 'hgnc_id', 'type', 'mode', 'mech',
             'syndrome']) + "\n"
        geneline = ("\t").join(
            ['4', '8846076', '8871839', 'HMX1', '5017', 'Probable DD gene',
             'Biallelic', 'Loss of function', 'OCULOAURICULAR SYNDROME']) + "\n"
        self.tempfile = tempfile.NamedTemporaryFile(mode="w")
        self.path = self.tempfile.name
        self.tempfile.write(header)
        self.tempfile.write(geneline)
        self.tempfile.flush()

    def test_load_genes(self):
        '''load genes from file'''
        self.assertEqual(load_genes(self.path), {
            '5017': {'chr': '4', 'start': '8846076', 'end': '8871839',
                     'symbol': 'HMX1', 'mode': {'Biallelic'},
                     'mechanism': {'Loss of function'},
                     'status': {'Probable DD gene'}}})

if __name__ == '__main__':
    unittest.main()