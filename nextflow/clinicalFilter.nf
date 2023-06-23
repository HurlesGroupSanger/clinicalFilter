nextflow.enable.dsl=2
include {CREATE_PED} from './modules/clinicalFilter.nf' 
include {CF} from './modules/clinicalFilter.nf' 

workflow {

	Channel.fromPath(params.probands)
	.splitCsv(header: false)
	.flatten()
	.map{row->tuple(row,
	 row.replaceAll("DDDP", "").substring(0,2),
	 row.replaceAll("DDDP", "").substring(2,4),
	 row.replaceAll("DDDP", "").substring(4,6))
	}
	.take(-1) // -1 means all samples
	.set{ch_probands}


	ch_ped = CREATE_PED(ch_probands, params.ped)

	ch_cf = CF(ch_ped, params.known_genes)
	ch_cf.view()
}