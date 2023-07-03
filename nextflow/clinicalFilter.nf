nextflow.enable.dsl=2
include {CREATE_PED} from './modules/clinicalFilter.nf' 
include {CF} from './modules/clinicalFilter.nf' 
include {CONCAT_RESULTS} from './modules/clinicalFilter.nf'

workflow {


// Starts from a list of individual proband PED files
if (params.containsKey("ped_list")) {

	Channel.fromPath(params.ped_list)
	.splitCsv(header: false, sep : "\t")
	.map{row->tuple(row[0],
	 row[0].replaceAll("DDDP", "").substring(0,2),
	 row[0].replaceAll("DDDP", "").substring(2,4),
	 row[0].replaceAll("DDDP", "").substring(4,6),
	 row[1])
	}
	.take(-1) // -1 means all samples
	.set{ch_ped}

}
else
// Starts from a proband list and a a study wide PED file
{
	// TODO : find a clever way not to have to run one process per proband (slow and demanding)
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

}

	ch_cf = CF(ch_ped, params.known_genes)

	ch_concat = CONCAT_RESULTS(ch_cf.map { tuple -> tuple[1] }.collect())
}
