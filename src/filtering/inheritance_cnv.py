import logging

from utils.utils import add_single_var_to_candidates
from utils.utils import add_compound_het_to_candidates

class CNVFiltering(object):

    def __init__(self, variants, family, genes, regions, trusted_variants, candidate_variants):
        self.variants = variants
        self.family = family
        self.genes = genes
        self.regions = regions
        self.trusted_variants = trusted_variants
        self.candidate_variants = candidate_variants
        self.parents = None
        if self.family.has_both_parents():
            self.parents = 'both'
        elif self.family.has_no_parents():
            self.parents = 'none'
        elif self.family.has_dad():
            self.parents = 'dad_only'
        elif self.family.has_mum():
            self.parents = 'mum_only'
        self.mum_aff = self.family.mum.affected
        self.dad_aff = self.family.dad.affected
        self.ddg2p_modes = set()

    def cnv_filter(self):
        # filter CNVs here
        if self.genes:
            self.cnv_filter_genes()
        if self.regions:
            # todo filtering for regions
            pass
        if self.trusted_variants:
            # todo filtering for trusted variants
            pass

    def cnv_filter_genes(self):
        if self.parents == 'both':
            self.cnv_filter_parents()
        elif self.parents == 'none':
            self.cnv_filter_no_parents()
        else:
            self.cnv_filter_single_parent()

    def cnv_filter_parents(self):
        for v in self.variants['child'].keys():
            if not self.variants['child'][v].is_cnv():
                continue
            if self.variants['child'][v].cnv_filter == 'Fail':
                logging.info(v + " CNV failed quality filters")
                continue
            #get gene modes
            modes = set()
            if self.genes:
                self.get_ddg2p_modes(v, modes)
        #inheritance
            self.inhmatch = self.cnv_inheritance_filter(v)
            if not self.inhmatch:
                #possible compound het
                posscomphet = self.cnv_candidate_compound_het_filter(v, modes)
                if not posscomphet:
                    logging.info(v + " failed CNV filter, inheritance doesn't "
                                     "match and not possible compound het")
            else:
                #non-ddg2p filter
                passnonddg2p = self.cnv_non_ddg2p_filter(v)
                if not passnonddg2p:
                    #ddg2p filter
                    passddg2p = self.cnv_ddg2p_filter(v)
                    if not passddg2p:
                        posscomphet = self.cnv_candidate_compound_het_filter(v, modes)
                        if not posscomphet:
                            logging.info(
                                v + " failed CNV filter, inheritance doesn't "
                                    "match and not possible compound het")

    def cnv_filter_no_parents(self):
        # non-ddg2p filter
        for v in self.variants['child'].keys():
            if not self.variants['child'][v].is_cnv():
                continue
            if self.variants['child'][v].cnv_filter == 'Fail':
                logging.info(v + " CNV failed quality filters")
                continue
            passnonddg2p = self.cnv_non_ddg2p_filter(v)
            if not passnonddg2p:
                # ddg2p filter
                passddg2p = self.cnv_ddg2p_filter(v)
                if not passddg2p:
                    posscomphet = self.cnv_candidate_compound_het_filter(v)
                    if not posscomphet:
                        logging.info(
                            v + " failed CNV filter, inheritance doesn't "
                                "match and not possible compound het")

    def cnv_filter_single_parent(self):
        pass

    def cnv_inheritance_filter(self, varid):
        # return True or False for pass or fail
        # inheritance matches parental affected status if:
        # paternal and father affected
        # OR
        # maternal and mum affected
        # OR
        # biparental plus CN = 0
        # OR
        # biparental and both/either parent(s) affected
        # OR
        # male (XY) proband X chromosome, maternal inh and mum unaffected
        # no variants are added to candidates at this stage
        if not self.variants['child'][varid].cnv_inh in ['paternal_inh', 'maternal_inh', 'biparental_inh']:
            return True
        if self.variants['child'][varid].cnv_inh == 'paternal_inh' and self.dad_aff:
            return True
        elif self.variants['child'][varid].cnv_inh == 'maternal_inh' and self.mum_aff:
            return True
        elif self.variants['child'][varid].cnv_inh == 'biparental_inh' and self.variants['child'][varid].cn == '0':
            return True
        elif self.variants['child'][varid].cnv_inh == 'biparental_inh' and (self.mum_aff or self.dad_aff):
            return True
        elif self.family.proband.sex == 'XY' and self.variants['child'][varid].chrom == 'X' and not self.mum_aff and self.variants['child'][varid].cnv_inh == 'maternal_inh':
            return True
        else:
            return False

    def cnv_candidate_compound_het_filter(self, varid, modes):
        # return True or False for pass or fail
        # could the CNV be part of a compound het? If so, add to candidate compound hets
        if not (self.variants['child'][varid].cn == '0' or self.variants['child'][varid].cn == '3'):
            return False
        else:
            # is any gene covered by the CNV biallelic, or CN = 1 and male hemizygous
            if 'Biallelic' in modes:
                add_compound_het_to_candidates(varid, self.variants['child'][varid], "-",
                                               "Biallelic",
                                               self.candidate_variants)
                return True
            elif 'Hemizygous' in modes and self.variants['child'][varid].cn == '1' and self.family.proband.sex == 'XY':
                add_compound_het_to_candidates(varid, self.variants['child'][varid], "-",
                                               "Hemizygous",
                                               self.candidate_variants)
                return True
            else:
                return False

    def cnv_non_ddg2p_filter(self, varid):
        # return True or False for pass or fail
        if int(self.variants['child'][varid].cnv_length) > 1000000:
            add_single_var_to_candidates(varid, self.variants['child'][varid], '-', '-',
                                         self.candidate_variants)
            return True
        else:
            return False

    def cnv_ddg2p_filter(self, varid):
        # return True or False for pass or fail
        cnvpass = False
        # get all genes covered by the CNV and go through each one at a time to see if any pass
        hgncids = self.variants['child'][varid].hgnc_id_all.split("|")
        for hid in hgncids:
            hgncid = hid[5:]
            if not hgncid in self.genes.keys():
                continue
            #fail duplications completely surrounding surround monoallelic, hemizygous and x-linked dominant genes with loss of function mechanism
            if self.variants['child'][varid].alt == "<DUP>" and self.variants['child'][varid].chrom == self.genes[hgncid]['chr']:
                #check modes
                dupmodes = set({"Monoallelic", "Hemizygous", "X-linked dominant"})
                if "Loss of function" in self.genes[hgncid]['mechanism'] and len(set.intersection(self.genes[hgncid]['mode'], dupmodes)) > 0:
                    if (int(self.variants['child'][varid].pos) < int(self.genes[hgncid]['start'])) and (int(self.variants['child'][varid].cnv_end) > int(self.genes[hgncid]['end'])):
                        logging.info(
                            varid + " duplication completely surrounds monoallelic, hemizygous or heterozygus gene with LoF mechanism")
            #Biallelic gene pass if copy number (CN) = 0 and mechanism in "Uncertain", "Loss of function", "Dominant negative"
            if int(self.variants['child'][varid].cn) >= 0 and "Biallelic" in self.genes[hgncid]['mode']:
                biallelicmechs = set({"Uncertain", "Loss of function", "Dominant negative"})
                if len(set.intersection(self.genes[hgncid]['mechanism'], biallelicmechs)) > 0:
                    cnvpass = True
                    add_single_var_to_candidates(varid,
                                                 self.variants['child'][varid],
                                                 hgncid, 'biallelic',
                                                 self.candidate_variants)
            #Monoallelic, X-linked dominant or Hemizygous in male pass if CN=0, 1 or 3 and any mechanism
            cns_wanted = ['0', '1', '3']
            if "Monoallelic" in self.genes[hgncid]['mode'] or "X-linked dominant" in self.genes[hgncid]['mode']:
                if self.variants['child'][varid].cn in cns_wanted:
                    cnvpass = True
                    add_single_var_to_candidates(varid,
                                                 self.variants['child'][varid],
                                                 hgncid, (",").join(self.genes[hgncid]['mode']),
                                                 self.candidate_variants)
            if "Hemizygous" in self.genes[hgncid]['mode'] and self.family.proband.sex == 'XY':
                cnvpass = True
                add_single_var_to_candidates(varid,
                                             self.variants['child'][varid],
                                             hgncid, (",").join(
                        self.genes[hgncid]['mode']),
                                             self.candidate_variants)
            #Hemizygous in female pass if CN=3 and mechanism = "Increased gene dosage"
            if "Hemizygous" in self.genes[
                hgncid]['mode'] and self.family.proband.sex == 'XX' and "Increased gene dosage" in self.genes[hgncid]['mechanism']:
                cnvpass = True
                add_single_var_to_candidates(varid,
                                             self.variants['child'][varid],
                                             hgncid, "Hemizygous",
                                             self.candidate_variants)
            #Pass intragenic DUP in monoallelic or X-linked dominant gene with loss of function mechanism and any part of the gene is outside of the CNV boundary
            if self.variants['child'][varid].alt == "<DUP>":
                if ("Monoallelic" in self.genes[hgncid]['mode'] or "X-linked dominant" in self.genes[hgncid]['mode']) and "Loss of function" in self.genes[hgncid]['mechanism']:
                    if int(self.variants['child'][varid].pos) > int(self.genes[hgncid]['start']) or int(self.variants['child'][varid].cnv_end) > int(self.genes[hgncid]['end']):
                        cnvpass = True
                        add_single_var_to_candidates(varid,
                                                     self.variants['child'][varid],
                                                     hgncid, (",").join(self.genes[hgncid]['mode']),
                                                     self.candidate_variants)

        return cnvpass

    def get_ddg2p_modes(self, varid, modes):
        # when a gene list is given, find the modes of all
        hgncids = self.variants['child'][varid].hgnc_id_all.split("|")
        for hid in hgncids:
            hgncid = hid[5:]
            if hgncid in self.genes.keys():
                genemodes = self.genes[hgncid]['mode']
                for m in genemodes:
                    modes.add(m)


