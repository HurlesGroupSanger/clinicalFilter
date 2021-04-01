"""copyright"""

# create inheritance report and populate for matrix creation for all
# combinations of parent and child gt, affected status and gene mode

def create_blank_inheritance_report():
    # create a blank report dict to report variants observed in all combinations
    # of trio genotype, affected status and gene mode
    inheritance_report = {}

    gt_matrix = {
        'dad_0/0_mum_0/0': 0,
        'dad_0/0_mum_0/1': 0,
        'dad_0/0_mum_1/1': 0,
        'dad_0/1_mum_0/0': 0,
        'dad_0/1_mum_0/1': 0,
        'dad_0/1_mum_1/1': 0,
        'dad_1/1_mum_0/0': 0,
        'dad_1/1_mum_0/1': 0,
        'dad_1/1_mum_1/1': 0
    }

    aff_matrix = {
        'dad_unaffected':{
            'mum_unaffected':gt_matrix.copy(),
            'mum_affected':gt_matrix.copy()
        },
        'dad_affected':{
            'mum_unaffected':gt_matrix.copy(),
            'mum_affected':gt_matrix.copy()
        }
    }

    inheritance_report['autosomal'] = {}
    inheritance_report['allosomal'] = {}

    inheritance_report['autosomal']['biallelic'] = {}
    inheritance_report['autosomal']['biallelic']['child_het'] = aff_matrix.copy()
    inheritance_report['autosomal']['biallelic']['child_hom'] = aff_matrix.copy()
    inheritance_report['autosomal']['monoallelic'] = {}
    inheritance_report['autosomal']['monoallelic']['child_het'] = aff_matrix.copy()
    inheritance_report['autosomal']['monoallelic']['child_hom'] = aff_matrix.copy()
    inheritance_report['autosomal']['mosaic'] = {}
    inheritance_report['autosomal']['mosaic']['child_het'] = aff_matrix.copy()
    inheritance_report['autosomal']['mosaic']['child_hom'] = aff_matrix.copy()
    inheritance_report['autosomal']['imprinted'] = {}
    inheritance_report['autosomal']['imprinted']['child_het'] = aff_matrix.copy()
    inheritance_report['autosomal']['imprinted']['child_hom'] = aff_matrix.copy()

    inheritance_report['allosomal']['hemizygous'] = {}
    inheritance_report['allosomal']['hemizygous']['child_het'] = aff_matrix.copy()
    inheritance_report['allosomal']['hemizygous']['child_hemi'] = aff_matrix.copy()
    inheritance_report['allosomal']['hemizygous']['child_hom'] = aff_matrix.copy()
    inheritance_report['allosomal']['X-linked_dominant'] = {}
    inheritance_report['allosomal']['X-linked_dominant']['child_het'] = aff_matrix.copy()
    inheritance_report['allosomal']['X-linked_dominant']['child_hemi'] = aff_matrix.copy()
    inheritance_report['allosomal']['X-linked_dominant']['child_hom'] = aff_matrix.copy()
    inheritance_report['allosomal']['X-linked_over_dominant'] = {}
    inheritance_report['allosomal']['X-linked_over_dominant']['child_het'] = aff_matrix.copy()
    inheritance_report['allosomal']['X-linked_over_dominant']['child_hemi'] = aff_matrix.copy()
    inheritance_report['allosomal']['X-linked_over_dominant']['child_hom'] = aff_matrix.copy()

    return inheritance_report

def populate_inheritance_report(inheritance_report, chromtype, mode, child_gt, mum_gt, dad_gt, mum_aff, dad_aff):
    '''populate the inheritance report - both parents present'''
    if chromtype == 'autosomal':
        if child_gt == '0/1':
            child_geno = 'child_het'
        elif child_gt == '1/1':
            child_geno = 'child_hom'
        else:
            print("Error - unrecocgised child GT")
            exit(1)

        if mum_aff == True:
            mum_aff_state = 'mum_affected'
        elif mum_aff == False:
            mum_aff_state = 'mum_unaffected'
        else:
            print("Error - invalid mum affected status")
            exit(1)

        if dad_aff == True:
            dad_aff_state = 'dad_affected'
        elif dad_aff == False:
            dad_aff_state = 'dad_unaffected'
        else:
            print("Error - invalid dad affected status")
            exit(1)

        parent_gts = "dad_" + dad_gt + "_mum_" + mum_gt
        inheritance_report[chromtype][mode][child_geno][dad_aff_state][mum_aff_state][parent_gts]+=1

    else:
        # todo allosomal and return error for unknown chromtype
        pass