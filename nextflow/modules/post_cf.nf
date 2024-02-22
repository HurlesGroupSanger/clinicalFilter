
process ANNOTATE_RESULTS {


	publishDir "${params.publish_dir}/${params.date}/", mode: "copy", pattern: "clinical_filter_results_*"

	input:
	path cf_results
	path id_mapping
	path previous_gene_list
	path decipher_variants_info
	path b37_cf_results
	path b38_cf_previous_results

	output:
	tuple path("clinical_filter_results_annotated.tsv"), path("clinical_filter_results_annotated_new_variants.tsv")

	script :
	"""
	cat > "post_cf.conf" <<EOF
	{
	"id_mapping": "$id_mapping",
	"previous_gene_list": "$previous_gene_list",
	"decipher_variants_info": "$decipher_variants_info",
	"b37_cf_results":"$b37_cf_results",
	"b38_cf_previous_results": "$b38_cf_previous_results",
	"latest_cf_results": "$cf_results",
	"outdir": "./"
	}
	EOF

	annotate_results.py post_cf.conf

	"""

}

