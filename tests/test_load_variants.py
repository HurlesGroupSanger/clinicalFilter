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

from variants.snv import SNV
from variants.cnv import CNV
from file_loading.load_vcfs import readvcf


class TestLoadVariants(unittest.TestCase):
    """make temporary VCF and test loading variants"""

    def setUp(self):
        ''' construct the bare minimum of lines for a VCF file to enable test
        load of variants'''
        self.maxDiff = None
        self.vcfheader = "##fileformat=VCFv4.2\n" + \
                    '##contig=<ID=1,length=248956422,assembly=GRCh38>' + "\n" + \
                    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">' + "\n" + \
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' + "\n" + \
                    '##FORMAT=<ID=AD,Number=1,Type=String,Description="AD">' + "\n" + \
                    '##FORMAT=<ID=PID,Number=1,Type=String,Description="PID">' + "\n" + \
                    '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="CN">' + "\n" + \
                    '##FORMAT=<ID=CIFER_INHERITANCE,Number=1,Type=String,Description="CIFER">' + "\n" + \
                    '##INFO=<ID=Consequence,Number=A,Type=String,Description="Consequence annotations from Ensembl VEP">' + "\n" + \
                    '##INFO=<ID=SYMBOL,Number=A,Type=String,Description="SYMBOL annotations from Ensembl VEP">' + "\n" + \
                    '##INFO=<ID=Gene,Number=A,Type=String,Description="Gene annotations from Ensembl VEP">' + "\n" + \
                    '##INFO=<ID=Feature,Number=A,Type=String,Description="Feature annotations from Ensembl VEP">' + "\n" + \
                    '##INFO=<ID=CANONICAL,Number=A,Type=String,Description="CANONICAL annotations from Ensembl VEP">' + "\n" + \
                    '##INFO=<ID=MANE,Number=A,Type=String,Description="MANE annotations from Ensembl VEP">' + "\n" + \
                    '##INFO=<ID=HGNC_ID,Number=A,Type=String,Description="HGNC_ID annotations from Ensembl VEP">' + "\n" + \
                    '##INFO=<ID=MAX_AF,Number=A,Type=Float,Description="MAX_AF annotations from Ensembl VEP">' + "\n" + \
                    '##INFO=<ID=MAX_AF_POPS,Number=A,Type=String,Description="MAX_AF_POPS annotations from Ensembl VEP">' + "\n" + \
                    '##INFO=<ID=DDD_AF,Number=A,Type=Float,Description="DDD unaffected parent allele frequency">' + "\n" + \
                    '##INFO=<ID=REVEL,Number=A,Type=Float,Description="Annotattion from REVEL">' + "\n" + \
                    '##INFO=<ID=PolyPhen,Number=A,Type=String,Description="PolyPhen annotations from Ensembl VEP">' + "\n" + \
                    '##INFO=<ID=HGVSc,Number=A,Type=String,Description="HGVSc annotations from Ensembl VEP">' + "\n" + \
                    '##INFO=<ID=HGVSp,Number=A,Type=String,Description="HGVSp annotations from Ensembl VEP">' + "\n" + \
                    '##INFO=<ID=CNVFILTER,Number=1,Type=String,Description="CNV FILTER">' + "\n" + \
                    '##INFO=<ID=HGNC_ID_ALL,Number=A,Type=String,Description="HGNC_ID_ALL annotations from Ensembl VEP">' + "\n" + \
                    '##INFO=<ID=SYMBOL_ALL,Number=A,Type=String,Description="SYMBOL_ALL annotations from Ensembl VEP">' + "\n" + \
                    '#CHROM' + "\t" + 'POS' + "\t" + 'ID' + "\t" + 'REF' + "\t" + 'ALT' + "\t" + 'QUAL' + "\t" + 'FILTER' + \
                    "\t" + 'INFO' + "\t" + 'FORMAT' + "\t" + 'sample1' + "\n"
        self.infofields = (";").join(
            ['Consequence=missense_variant', 'SYMBOL=MECP1', 'Gene=ENSG01234',
             'Feature=ENST01234', 'CANONICAL=YES', 'MANE=NM01234',
             'HGNC_ID=HGNC:123', 'MAX_AF=0.0003', 'MAX_AF_POPS=.',
             'DDD_AF=0.0001',
             'REVEL=0.8', 'PolyPhen=Benign', 'HGVSc=string', 'HGVSp=string'])
        self.variantline = ("\t").join(
            ['1', '1339911', 'rs1234', 'A', 'G', '.', '.', self.infofields, 'GT:GQ:AD:PID',
             '1/1:99:.:.']) + "\n"
        self.commoninfofields = (";").join(
            ['Consequence=missense_variant', 'SYMBOL=MECP1', 'Gene=ENSG01234',
             'Feature=ENST01234', 'CANONICAL=YES', 'MANE=NM01234',
             'HGNC_ID=HGNC:123', 'MAX_AF=0.3', 'MAX_AF_POPS=.',
             'DDD_AF=0.0001',
             'REVEL=0.8', 'PolyPhen=Benign', 'HGVSc=string', 'HGVSp=string'])
        self.commonvariantline = ("\t").join(
            ['1', '1339915', 'rs1234', 'A', 'G', '.', '.', self.commoninfofields,
             'GT:GQ:AD:PID',
             '1/1:99:.:.']) + "\n"
        self.var3infofields = (";").join(
            ['Consequence=missense_variant', 'SYMBOL=MECP1', 'Gene=ENSG01234',
             'Feature=ENST01234', 'CANONICAL=YES', 'MANE=NM01234',
             'HGNC_ID=HGNC:123', 'MAX_AF=0.0003', 'MAX_AF_POPS=.',
             'DDD_AF=0.0001', 'REVEL=0.8', 'PolyPhen=Benign', 'HGVSc=string',
             'HGVSp=string'])
        self.var3variantline = ("\t").join(
            ['1', '1449915', 'rs1234', 'A', 'G', '.', '.', self.var3infofields,
             'GT:GQ:AD:PID',
             '1/1:99:.:.']) + "\n"
        self.cnvinfofields = (";").join(
            ['Consequence=transcript_ablation', 'SYMBOL=MECP1', 'Gene=ENSG01234',
             'Feature=ENST01234', 'CANONICAL=YES', 'MANE=NM01234',
             'HGNC_ID=HGNC:123', 'CNVFILTER=Pass', 'SYMBOL_ALL=MECP1|MECP2|MECP3',
             'HGNC_ID_ALL=HGNC:1|HGNC:2|HGNC:3'])
        self.cnvline = ("\t").join(['1', '123456', '.', 'T', '<DEL>', '.', '.',
                        self.cnvinfofields, 'CN:CIFER_INHERITANCE', '1:maternal_inh'])

    def test_load_snv(self):
        '''load two variant lines - the common variant is ignored'''
        self.tempfile = tempfile.NamedTemporaryFile(mode="w")
        self.path = self.tempfile.name
        self.tempfile.write(self.vcfheader)
        self.tempfile.write(self.variantline)
        self.tempfile.write(self.commonvariantline)
        self.tempfile.write(self.var3variantline)
        self.tempfile.flush()

        self.assertEqual(readvcf(self.path, None, 'XY'), {'1_1339911_A_G': SNV(
            {'chrom': '1', 'pos': '1339911', 'ref': 'A', 'alt': 'G',
             'consequence': 'missense_variant', 'ensg': 'ENSG01234',
             'symbol': 'MECP1', 'feature': 'ENST01234', 'canonical': 'YES',
             'mane': 'NM01234', 'hgnc_id': 'HGNC:123', 'max_af': '0.0003',
             'max_af_pops': '.', 'ddd_af': '0.0001', 'revel': '0.8',
             'polyphen': 'Benign', 'hgvsc': 'string', 'hgvsp': 'string',
             'dnm': '.', 'gt': '1/1',
             'gq': '99', 'sex': 'XY'}), '1_1449915_A_G': SNV(
            {'chrom': '1', 'pos': '1449915', 'ref': 'A',
             'alt': 'G', 'consequence': 'missense_variant',
             'ensg': 'ENSG01234', 'symbol': 'MECP1',
             'feature': 'ENST01234', 'canonical': 'YES',
             'mane': 'NM01234', 'hgnc_id': 'HGNC:123',
             'max_af': '0.0003', 'max_af_pops': '.',
             'ddd_af': '0.0001', 'revel': '0.8',
             'polyphen': 'Benign', 'hgvsc': 'string',
             'hgvsp': 'string', 'dnm':'.', 'gt': '1/1', 'gq': '99',
             'sex': 'XY'})})

    def test_load_cnv(self):
        '''load CNV'''
        self.tempfile = tempfile.NamedTemporaryFile(mode="w")
        self.path = self.tempfile.name
        self.tempfile.write(self.vcfheader)
        self.tempfile.write(self.cnvline)
        self.tempfile.flush()

        self.assertEqual(readvcf(self.path, None, 'XY'), {'1_123456_T_<DEL>': CNV(
            {'chrom': '1', 'pos': '123456', 'ref': 'T', 'alt': '<DEL>',
             'consequence': 'transcript_ablation', 'ensg': 'ENSG01234',
             'symbol': 'MECP1', 'feature': 'ENST01234', 'canonical': 'YES',
             'mane': 'NM01234', 'hgnc_id': 'HGNC:123', 'cnv_filter': 'Pass',
             'hgcn_id_all': 'HGNC:1|HGNC:2|HGNC:3',
             'symbol_all': 'MECP1|MECP2|MECP3', 'sex': 'XY', 'cn': '1',
             'cnv_inh':'maternal_inh'})})

if __name__ == '__main__':
    unittest.main()
