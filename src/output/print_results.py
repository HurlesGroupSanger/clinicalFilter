"""copyright"""

# -identify variants in both compound het and single variants dicts
# -identify and flag possible MNVs
# -flag variants in cis in monoallelic genes

def create_output(candidate_variants, family, outdir):
    identify_mnvs(candidate_variants)
    identify_variants_in_cis(candidate_variants)
    print_output(candidate_variants, family, outdir)

def identify_mnvs(candidate_variants):
    pass
# check both phasing and parental information 

def identify_variants_in_cis(candidate_variants):
    pass

def print_output(candidate_variants, family, outdir):
    pass