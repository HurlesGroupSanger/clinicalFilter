# Installation

git clone https://github.com/RuthEberhardt/clinicalFilter.git

Add the src directory from the cloned repository to your pythonpath, **DIR** represents the dirctory the code is cloned into:

export PYTHONPATH=**DIR**/clinicalFilter/src/:$PYTHONPATH

# Running clinicai filtering

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
