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


class Variant(object):
    """
    Generic variant class, inherited by more specific classes such as CNV
     and SNV
     """

    def __init__(self, vardata):
        for key in vardata:
            setattr(self, key, vardata[key])

        self.genotype = None
        self.triogenotype = None
        self.set_genotype()
        self.standardise_chromosome()
        self.parse_hgnc_id()

    def __eq__(self, other):
        return self.chrom == other.chrom and \
               self.pos == other.pos and \
               self.ref == other.ref and \
               self.alt == other.alt

    def standardise_chromosome(self):
        """
        Ensure chromosome is 1-22,X,Y
        """
        if self.chrom.startswith('Chr') or self.chrom.startswith('chr'):
            self.chrom = self.chrom[3:]

    def parse_hgnc_id(self):
        """
        Strip HGNC: from hgnc_id
        """
        if self.hgnc_id.startswith('HGNC:'):
            self.hgnc_id = self.hgnc_id[5:]





