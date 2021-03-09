"""
copyright
"""

from file_loading.load_genes_and_regions import load_genes
from file_loading.loadvcfs import load_variants
from variants.triogenotype import add_trio_genotypes
from filtering.preinheritance_filtering import preinheritance_filter
from filtering.inheritance_filtering import inheritance_filter
from filtering.postinheritance_filter import postinheritance_filter

def filter_trio(family, genes_file, regions_file, trusted_variants_file, outdir):
    """filter each trio"""

    # if genes, regions or variants files are present we can create a list of
    # regions to load and therefore load fewer variants
    vcfregions = set()
    genes = None
    regions = None
    trusted_variants = None

    if genes_file:
        genes = load_genes(genes_file)
        for g in genes.keys():
            reg = genes[g]['chr'] + "\t" + genes[g]['start'] + "\t" + genes[g]['end']
            vcfregions.add(reg)

    if regions_file:
        #TODO add regions to the vcfregions set and populate regions variable
        pass

    if trusted_variants_file:
        #TODO add trusted variants locations to the vcfregions set and populate
        # trusted_regions variable
        pass

    if len(vcfregions) > 0:
        variants = load_variants(family, outdir, vcfregions)
    else:
        variants = load_variants(family, outdir)

    #add trio genotypes for each variant
    add_trio_genotypes(family, variants)

    #preinheritance filters
    variants_per_gene = preinheritance_filter(variants)
    # print(variants_per_gene)

    #inheritance filters
    candidate_variants = inheritance_filter(variants_per_gene, family, genes, regions, trusted_variants)
    # print(candidate_variants['compound_hets'])
    # print(candidate_variants['single_variants'])

    #post inheritance filters
    postinheritance_filter(candidate_variants, family)
    print(candidate_variants['compound_hets'].keys())
    print(candidate_variants['single_variants'].keys())

