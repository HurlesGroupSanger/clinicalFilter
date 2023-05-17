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

import unittest
import tempfile

from file_loading.load_genes_and_regions import load_genes


class TestLoadGenesRegions(unittest.TestCase):
    def setUp(self):
        """create genes file to load"""
        self.maxDiff = None
        header = ("\t").join(["chr", "start", "stop", "gene", "hgnc_id", "type", "mode", "mech", "syndrome"]) + "\n"
        geneline = ("\t").join(
            [
                "4",
                "8846076",
                "8871839",
                "HMX1",
                "5017",
                "Probable DD gene",
                "Biallelic",
                "Loss of function",
                "OCULOAURICULAR SYNDROME",
            ]
        ) + "\n"
        self.tempfile = tempfile.NamedTemporaryFile(mode="w")
        self.path = self.tempfile.name
        self.tempfile.write(header)
        self.tempfile.write(geneline)
        self.tempfile.flush()

    def test_load_genes(self):
        """load genes from file"""
        self.assertEqual(
            load_genes(self.path),
            {
                "5017": {
                    "chr": "4",
                    "start": "8846076",
                    "end": "8871839",
                    "symbol": "HMX1",
                    "mode": {"Biallelic"},
                    "mechanism": {"Loss of function"},
                    "status": {"Probable DD gene"},
                }
            },
        )


if __name__ == "__main__":
    unittest.main()
