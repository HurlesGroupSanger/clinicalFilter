

process CREATE_PED {
	

	publishDir "${params.publish_dir}/${params.date}/", mode: 'copy', pattern : "list_ped.tsv"


	input:
    path(list_vcfs)
    path(list_probands)
	path(families_ped)
	val(publish_dir)
	val(date)

	output:
	path("list_ped.tsv")
	
	script :
	"""
	create_ped.py $list_vcfs $list_probands $families_ped $publish_dir/$date 
	"""
}

process CF {
	
    tag "CF_$stable_id"

	publishDir "${params.publish_dir}/${params.date}/DD/DP/$ab/$cd/$ef/${stable_id}/clinical_filter/", mode: 'copy', pattern: "${stable_id}_clinical_filter_inheritance_report.txt"
	publishDir "${params.publish_dir}/${params.date}/DD/DP/$ab/$cd/$ef/${stable_id}/clinical_filter/", mode: 'copy', pattern: "${stable_id}_clinical_filter.tsv"
	publishDir "${params.publish_dir}/${params.date}/DD/DP/$ab/$cd/$ef/${stable_id}/clinical_filter/", mode: 'copy', pattern: "${stable_id}_clinical_filter.log"


	beforeScript = params.useModules
		? """
		module load ${params.bcftoolsModule}
		export PYTHONPATH=${baseDir}/../src:\$PYTHONPATH
		"""
		: """
		export PYTHONPATH=${baseDir}/../src:\$PYTHONPATH
		"""

	input:
    tuple val(stable_id),val(ab),val(cd),val(ef), path(ped)
	path(known_genes)

	output :
	tuple path("${stable_id}_clinical_filter_inheritance_report.txt"), path("${stable_id}_clinical_filter.tsv"), path("${stable_id}_clinical_filter.log")

	script :
	"""
	runclinicalfiltering.py --ped $ped --known-genes $known_genes --outdir ./
	"""
}


process CONCAT_RESULTS {

	tag "CONCAT_RESULTS"

	publishDir "${params.publish_dir}/${params.date}/", mode: "copy", pattern: "clinical_filter_results.tsv"

	input:
	path clinical_filter_individual_results

	output:
	path "clinical_filter_results.tsv"


	script:
	"""
	awk 'FNR==1 && NR!=1{next;}{print}' $clinical_filter_individual_results > clinical_filter_results.tsv
	"""

}