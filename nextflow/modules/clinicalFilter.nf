

process CREATE_PED {
	
    tag "CREATE_PED_$stable_id"

	publishDir "${params.publish_dir}/DD/DP/$ab/$cd/$ef/${stable_id}/clinical_filter/", mode: 'copy', pattern: "${stable_id}.ped"


	input:
    tuple val(stable_id),val(ab),val(cd),val(ef)
    path(all_ped)

	output:
	tuple val(stable_id),val(ab),val(cd),val(ef), path("${stable_id}.ped")

	script :
	"""
	#!/usr/bin/env python
	import pandas as pd

	all_ped_df = pd.read_csv("${all_ped}", sep="\\t")
	family_id = all_ped_df.loc[all_ped_df.individual_id == "${stable_id}", "family_id"]
	assert(len(family_id) == 1)

	all_ped_df.loc[all_ped_df.family_id == family_id.iloc[0]].to_csv("${stable_id}.ped", sep="\\t", header=False, index=False)
	"""
}

process CF {
	
    tag "CF_$stable_id"

	publishDir "${params.publish_dir}/DD/DP/$ab/$cd/$ef/${stable_id}/clinical_filter/", mode: 'copy', pattern: "${stable_id}_clinical_filter_inheritance_report.txt"
	publishDir "${params.publish_dir}/DD/DP/$ab/$cd/$ef/${stable_id}/clinical_filter/", mode: 'copy', pattern: "${stable_id}_clinical_filter.txt"
	publishDir "${params.publish_dir}/DD/DP/$ab/$cd/$ef/${stable_id}/clinical_filter/", mode: 'copy', pattern: "${stable_id}_clinical_filter.log"

	beforeScript "export PYTHONPATH=${baseDir}/../src/:\$PYTHONPATH"

	input:
    tuple val(stable_id),val(ab),val(cd),val(ef), path(ped)
	path(known_genes)

	output :
	tuple path("${stable_id}_clinical_filter_inheritance_report.txt"), path("${stable_id}_clinical_filter.txt"), path("${stable_id}_clinical_filter.log")

	script :
	"""
	runclinicalfiltering.py --ped $ped --known-genes $known_genes --outdir ./
	"""
}


process CONCAT_RESULTS {

	tag "CONCAT_RESULTS"

	publishDir "${params.publish_dir}/", mode: "copy", pattern: "clinical_filter_results.tsv"

	input:
	path clinical_filter_individual_results

	output:
	path "clinical_filter_results.tsv"


	script:
	"""
	awk 'FNR==1 && NR!=1{next;}{print}' $clinical_filter_individual_results > clinical_filter_results.tsv
	"""

}