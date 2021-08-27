"""
copyright
"""

from file_loading.load_genes_and_regions import load_genes
from file_loading.load_vcfs import load_variants
from variants.trio_genotype import add_trio_genotypes
from filtering.preinheritance_filtering import PreInheritanceFiltering
from filtering.inheritance_filtering import InheritanceFiltering
from filtering.inheritance_cnv import CNVFiltering
from filtering.postinheritance_filter import PostInheritanceFiltering
from filtering.inheritance_report import InheritanceReport
from filtering.compound_hets import CompoundHetScreen

class Filter(object):
    '''Class for filtering variants'''

    def __init__(self, family, known_genes, known_regions,
                    trusted_variants, outdir):
        self.family = family
        self.known_genes = known_genes
        self.known_regions = known_regions
        self.trusted_variants = trusted_variants
        self.outdir = outdir
        self.candidate_variants = None
        self.candidate_variants = {}
        self.candidate_variants['single_variants'] = {}
        self.candidate_variants['compound_hets'] = {}
        self.inhreport = None
        self.inhreport = InheritanceReport()
        self.screened_candidate_variants = {}

    def filter_trio(self):
        """filter each trio"""

        # if genes, regions or variants files are present we can create a list of
        # regions to load and therefore load fewer variants
        vcfregions = set()
        genes = None
        regions = None
        trusted_variants = None

        if self.known_genes:
            genes = load_genes(self.known_genes)
            # for g in genes.keys():
            #     reg = genes[g]['chr'] + "\t" + genes[g]['start'] + "\t" + genes[g]['end']
            #     vcfregions.add(reg)

        if self.known_regions:
            #TODO add regions to the vcfregions set and populate regions variable
            pass

        if self.trusted_variants:
            #TODO add trusted variants locations to the vcfregions set and populate
            # trusted_regions variable
            pass

        if len(vcfregions) > 0:
            variants = load_variants(self.family, self.outdir, vcfregions)
        else:
            variants = load_variants(self.family, self.outdir)

        #add trio genotypes for each variant
        add_trio_genotypes(self.family, variants)

        #preinheritance filters
        preinheritancefilter = PreInheritanceFiltering(variants)
        variants_per_gene = preinheritancefilter.preinheritance_filter()

        #inheritance filters for SNVs
        inheritancefilter = InheritanceFiltering(variants_per_gene, self.family, genes, regions, trusted_variants, self.candidate_variants, self.inhreport)
        # candidate_variants, inheritance_report = inheritancefilter.inheritance_filter()
        inheritancefilter.inheritance_filter()

        # inheritance filters for CNVs
        cnvfilter = CNVFiltering(variants, self.family, genes, regions, trusted_variants, self.candidate_variants)
        cnvfilter.cnv_filter()

        # compound het screen
        compoundhets = CompoundHetScreen(self.candidate_variants, self.family)
        compoundhets.screen_compound_hets()
        # self.screened_candidate_variants = compoundhets.screen_compound_hets()

        #post inheritance filters
        postinheritancefilter = PostInheritanceFiltering(self.candidate_variants, self.family)
        filtered_candidate_variants = postinheritancefilter.postinheritance_filter()

        # print(candidate_variants)
        # print(filtered_candidate_variants)
        # print(candidate_variants['compound_hets'].keys())
        # print(candidate_variants['single_variants'].keys())

        return filtered_candidate_variants, self.inhreport

