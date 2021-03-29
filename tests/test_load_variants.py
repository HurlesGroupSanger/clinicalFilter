import unittest
import tempfile

from variants.snv import SNV
from file_loading.loadvcfs import readvcf


class TestLoadVariants(unittest.TestCase):
    """make temporary VCF and test loading variants"""

    def setUp(self):
        ''' construct the bare minimum of lines for a VCF file to enable test
        load of variants'''
        self.maxDiff = None
        vcfheader = "##fileformat=VCFv4.2\n" + \
                    '##contig=<ID=1,length=248956422,assembly=GRCh38>' + "\n" + \
                    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">' + "\n" + \
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' + "\n" + \
                    '##FORMAT=<ID=AD,Number=1,Type=String,Description="AD">' + "\n" + \
                    '##FORMAT=<ID=PID,Number=1,Type=String,Description="PID">' + "\n" + \
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
                    '#CHROM' + "\t" + 'POS' + "\t" + 'ID' + "\t" + 'REF' + "\t" + 'ALT' + "\t" + 'QUAL' + "\t" + 'FILTER' + \
                    "\t" + 'INFO' + "\t" + 'FORMAT' + "\t" + 'sample1' + "\n"
        infofields = (";").join(
            ['Consequence=missense_variant', 'SYMBOL=MECP1', 'Gene=ENSG01234',
             'Feature=ENST01234', 'CANONICAL=YES', 'MANE=NM01234',
             'HGNC_ID=HGNC:123', 'MAX_AF=0.0003', 'MAX_AF_POPS=.',
             'DDD_AF=0.0001',
             'REVEL=0.8', 'PolyPhen=Benign', 'HGVSc=string', 'HGVSp=string'])
        variantline = ("\t").join(
            ['1', '1339911', 'rs1234', 'A', 'G', '.', '.', infofields, 'GT:GQ:AD:PID',
             '1/1:99:.:.']) + "\n"
        commoninfofields = (";").join(
            ['Consequence=missense_variant', 'SYMBOL=MECP1', 'Gene=ENSG01234',
             'Feature=ENST01234', 'CANONICAL=YES', 'MANE=NM01234',
             'HGNC_ID=HGNC:123', 'MAX_AF=0.3', 'MAX_AF_POPS=.',
             'DDD_AF=0.0001',
             'REVEL=0.8', 'PolyPhen=Benign', 'HGVSc=string', 'HGVSp=string'])
        commonvariantline = ("\t").join(
            ['1', '1339915', 'rs1234', 'A', 'G', '.', '.', commoninfofields,
             'GT:GQ:AD:PID',
             '1/1:99:.:.']) + "\n"
        var3infofields = (";").join(
            ['Consequence=missense_variant', 'SYMBOL=MECP1', 'Gene=ENSG01234',
             'Feature=ENST01234', 'CANONICAL=YES', 'MANE=NM01234',
             'HGNC_ID=HGNC:123', 'MAX_AF=0.0003', 'MAX_AF_POPS=.',
             'DDD_AF=0.0001', 'REVEL=0.8', 'PolyPhen=Benign', 'HGVSc=string',
             'HGVSp=string'])
        var3variantline = ("\t").join(
            ['1', '1449915', 'rs1234', 'A', 'G', '.', '.', var3infofields,
             'GT:GQ:AD:PID',
             '1/1:99:.:.']) + "\n"
        self.tempfile = tempfile.NamedTemporaryFile(mode="w")
        self.path = self.tempfile.name
        self.tempfile.write(vcfheader)
        self.tempfile.write(variantline)
        self.tempfile.write(commonvariantline)
        self.tempfile.write(var3variantline)
        self.tempfile.flush()

        with open('/nfs/team29/re3/new_clinical_filtering/tmp.vcf', 'w') as o:
            o.write(vcfheader)
            o.write(variantline)
            o.write(commonvariantline)
            o.write(var3variantline)

    def test_load_variants(self):
        '''load two variant lines - the common variant is ignored'''
        self.assertEqual(readvcf(self.path, None, 'XY'), {'1_1339911_A_G': SNV(
            {'chrom': '1', 'pos': '1339911', 'ref': 'A', 'alt': 'G',
             'consequence': 'missense_variant', 'ensg': 'ENSG01234',
             'symbol': 'MECP1', 'feature': 'ENST01234', 'canonical': 'YES',
             'mane': 'NM01234', 'hgnc_id': 'HGNC:123', 'max_af': '0.0003',
             'max_af_pops': '.', 'ddd_af': '0.0001', 'revel': '0.8',
             'polyphen': 'Benign', 'hgvsc': 'string', 'hgvsp': 'string',
             'denovo_snv': False, 'denovo_indel': False, 'gt': '1/1',
             'gq': '99', 'sex': 'XY'}), '1_1449915_A_G': SNV(
            {'chrom': '1', 'pos': '1449915', 'ref': 'A',
             'alt': 'G', 'consequence': 'missense_variant',
             'ensg': 'ENSG01234', 'symbol': 'MECP1',
             'feature': 'ENST01234', 'canonical': 'YES',
             'mane': 'NM01234', 'hgnc_id': 'HGNC:123',
             'max_af': '0.0003', 'max_af_pops': '.',
             'ddd_af': '0.0001', 'revel': '0.8',
             'polyphen': 'Benign', 'hgvsc': 'string',
             'hgvsp': 'string', 'denovo_snv': False,
             'denovo_indel': False, 'gt': '1/1', 'gq': '99',
             'sex': 'XY'})})

if __name__ == '__main__':
    unittest.main()
