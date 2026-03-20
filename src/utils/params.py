# Decipher genotype values
DECIPHER_GENOTYPE_HOMOZYGOUS = "homozygous"
DECIPHER_GENOTYPE_HETEROZYGOUS = "heterozygous"
DECIPHER_GENOTYPE_HEMIZYGOUS = "hemizygous"
DECIPHER_GENOTYPE_NA = "."


# Decipher inheritance values
DECIPHER_INHERITANCE_BIPARENTAL = "biparental"
DECIPHER_INHERITANCE_DENOVO = "de_novo_confirmed"
DECIPHER_INHERITANCE_MATERNAL = "maternal"
DECIPHER_INHERITANCE_PATERNAL = "paternal"
DECIPHER_INHERITANCE_DENOVO_MOSAIC = "de_novo_mosaic"
DECIPHER_INHERITANCE_MATERNAL_MOSAIC = "maternal_mosaic"
DECIPHER_INHERITANCE_PATERNAL_MOSAIC = "paternal_mosaic"
DECIPHER_INHERITANCE_UNKNOWN = "unknown"
DECIPHER_INHERITANCE_NA = "."

# Threshold used on allelic depth to categorise the variant as mosaic for DECIPHER inheritance
THRESHOLD_AD_MOSAICITY = 0.3

# CEP thresholds
SPLICE_AI_THRESHOLD = 0.2
REVEL_THRESHOLD = 0.4
REVEL_THRESHOLD_COMPOUNDHET_SINGLETON = 0.7

# Functional consequences
MODERATE_HIGH_IMPACT_CONSEQUENCES = [
    "frameshift_variant",
    "missense_variant",
    "splice_donor_variant",
    "splice_acceptor_variant",
    "start_lost",
    "stop_gained",
    "protein_altering_variant",
    "transcript_ablation",
    "transcript_amplification",
    "inframe_insertion",
    "inframe_deletion",
    "stop_lost",
]

LOF_CONSEQUENCES = [
    "frameshift_variant",
    "splice_donor_variant",
    "splice_acceptor_variant",
    "stop_gained",
    "transcript_ablation",
]
