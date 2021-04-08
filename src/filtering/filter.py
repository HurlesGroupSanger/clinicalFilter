"""
copyright
"""

from file_loading.load_genes_and_regions import load_genes
from file_loading.loadvcfs import load_variants
from variants.triogenotype import add_trio_genotypes
from filtering.preinheritance_filtering import PreInheritanceFiltering
from filtering.inheritance_filtering import InheritanceFiltering
from filtering.postinheritance_filter import PostInheritanceFiltering

class Filter(object):
    '''Class for filtering variants'''

    def __init__(self, family, known_genes, known_regions,
                    trusted_variants, outdir):
        self.family = family
        self.known_genes = known_genes
        self.known_regions = known_regions
        self.trusted_variants = trusted_variants
        self.outdir = outdir

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
            for g in genes.keys():
                reg = genes[g]['chr'] + "\t" + genes[g]['start'] + "\t" + genes[g]['end']
                vcfregions.add(reg)

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
        Preinheritancefilter = PreInheritanceFiltering(variants)
        variants_per_gene = Preinheritancefilter.preinheritance_filter()

        #inheritance filters
        Inheritancefilter = InheritanceFiltering(variants_per_gene, self.family, genes, regions, trusted_variants)
        candidate_variants, inheritance_report = Inheritancefilter.inheritance_filter()
        # candidate_variants, inheritance_report = inheritance_filter(variants_per_gene, self.family, genes, regions, trusted_variants)
        # import pprint as pp
        # pp.pprint(inheritance_report)
        # print(candidate_variants)
        # exit(0)

        #post inheritance filters
        Postinheritancefilter = PostInheritanceFiltering(candidate_variants, self.family)
        filtered_candidate_variants = Postinheritancefilter.postinheritance_filter()

        # print(candidate_variants)
        # print(filtered_candidate_variants)
        # print(candidate_variants['compound_hets'].keys())
        # print(candidate_variants['single_variants'].keys())

        return filtered_candidate_variants, inheritance_report

