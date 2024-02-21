nextflow.enable.dsl=2
include {CREATE_PED} from './modules/clinicalFilter.nf' 
include {CF} from './modules/clinicalFilter.nf' 
include {CONCAT_RESULTS} from './modules/clinicalFilter.nf'

workflow {


ch_pedlist = CREATE_PED(params.list_vcfs, params.list_probands, params.families_ped, params.publish_dir, params.date)


ch_pedlist
	.splitCsv(header: false, sep : "\t")
	.map{row->tuple(row[0],
		row[0].replaceAll("DDDP", "").substring(0,2),
		row[0].replaceAll("DDDP", "").substring(2,4),
		row[0].replaceAll("DDDP", "").substring(4,6),
		row[1])
	}
	.take(-1) // -1 means all samples
	.set{ch_ped}


ch_cf = CF(ch_ped, params.known_genes)
ch_concat = CONCAT_RESULTS(ch_cf.map { tuple -> tuple[1] }.collect())

}
