# Installation

git clone https://github.com/RuthEberhardt/clinicalFilter.git

Add the src directory from the cloned repository to your pythonpath, **DIR** represents the dirctory the code is cloned into:

export PYTHONPATH=**DIR**/clinicalFilter/src/:$PYTHONPATH

# Requirements

Requires bctfools to be on the PATH to run

# Running clinical filtering

DIR represents the dirctory the code is cloned into

To view help
```sh
python3 DIR/clinicalFilter/runclinicalfiltering.py  --help
```

To run clinical filtering using a gene list from a ped file
```sh
python3 DIR/clinicalFilter/runclinicalfiltering.py \
--ped PED_PATH \
--known-genes GENES_FILE \
--outdir OUTPUT_DIR 
```

To run clinical filtering using a gene list without a ped file
```sh
python3 DIR/clinicalFilter/runclinicalfiltering.py \
--child PATH_TO_CHILD_VCF \
--sex CHILD SEX (ed XX,XY) \
--mother PATH_TO_MUM_VCF \
--father PATH_TO_DAD_VCF \
--mum_aff MUM_AFFECTED_STATUS (1 = unaffected, 2 = affected) \
--daf_aff DAD_AFFECTED_STATUS (1 = unaffected, 2 = affected) \
--known-genes GENES_FILE \
--outdir OUTPUT_DIR 
```

# Input files

**VCF files**

VCF files should be annotated with VEP (version 104 tested) including the REVEL plugin, and the VEP annotation split using bcftools split-vep plugin:
/software/ddd/external/bcftools/bcftools +split-vep -c - -s worst INPUT_VCF | bgzip -c > OUTPUT_VCF

The following VEP annotation is used:
Consequence, Gene, SYMBOL, Feature, CANONICAL, MANE_SELECT, MANE_PLUS_CLINICAL, HGNC_ID, MAX_AF, MAX_AF_POPS, REVEL, PolyPhen, Protein_position, HGVSc, HGVSp

Allele count annotation from gnomAD is required:
AC_XX, AN_XX, nhomalt_XX AC_XY AN_XY, nhomalt_XY

For trios, DNM annotation from the bcftools trio-dnm2 plugin, http://samtools.github.io/bcftools/howtos/plugin.trio-dnm2.html, is required. The following fields are used:
pp_trio_DNM2, pp_DNG

CNVs can be added from any caller if desired. If CNVs are present the following annotation should be present:
INFO: CNVFILTER (Pass or Fail)
INFO: END
FORMAT: CIFER_INHERITANCE (CNV inheritance: biparental_inh / maternal_inh / paternal_inh / not_inherited / unceetain / unable_to_evaluate_probes / false_positive)
FORMAT: CN (Copy number)

Optionally, cohort specific allele frequences can be added to the VCF files. Currently these should be added as DDD_AF (unaffected paretnal allele frquency) and DDD_father_AF (unaffected father allele frequency for X chromosome). If the cohort is too small these may cause all variants to fail allele frequncy thresholds.

**ped file**

Tab separated file with the following fields:
family id
person id
father id (0 is used when there is no father)
mother id (0 is used where there is no mother)
chromosomal sex, abnormal karyotypes are supported (eg XXY)
affected status (2 = affected, 1 = unaffected)
path to VCF for individual

The proband-list option can be used to load a subset of probands from a ped file

Example:
```sh
fam1	proband_id	dad_id	mum_id	XX	2	PATH_TO_PROBAND_VCF
fam1	mum_id	0	0	XX	1	PATH_TO_MUM_VCF
fam1	dad_id	0	0	XY	1	PATH_TO_DAD_VCF
```

**gene list** 

A tab-separated gene list in the following format:

```
chr	start	stop	gene	hgnc_id	type	mode	mech	syndrome
12	88049016	88142099	CEP290	29021	Confirmed DD gene	Biallelic	Loss of function	JOUBERT SYNDROME TYPE 5
15	74179466	74212267	STRA6	30650	Confirmed DD gene	Biallelic	Loss of function	MICROPHTHALMIA SYNDROMIC TYPE 9
15	74890005	74902219	MPI	7216	Confirmed DD gene	Biallelic	Loss of function	CONGENITAL DISORDERS OF GLYCOSYLATION
```
