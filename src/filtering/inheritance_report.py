"""
Copyright (c) 2021 Genome Research Limited
Author: Ruth Eberhardt <re3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import copy


class InheritanceReport(object):
    """
    Create inheritance report and populate for matrix creation for all
    combinations of parent and child gt, affected status and gene mode
    """

    def __init__(self):
        self.inheritance_report = self.create_blank_inheritance_report()

    def create_blank_inheritance_report(self):
        """
        Create a blank report dict to report variants observed in all
        combinations of trio genotype, affected status and gene mode
        """
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
            'dad_unaffected': {
                'mum_unaffected': copy.deepcopy(gt_matrix),
                'mum_affected': copy.deepcopy(gt_matrix)
            },
            'dad_affected': {
                'mum_unaffected': copy.deepcopy(gt_matrix),
                'mum_affected': copy.deepcopy(gt_matrix)
            }
        }

        inheritance_report['autosomal'] = {}
        inheritance_report['allosomal'] = {}

        inheritance_report['autosomal']['biallelic'] = {}
        inheritance_report['autosomal']['biallelic'][
            'child_het'] = copy.deepcopy(aff_matrix)
        inheritance_report['autosomal']['biallelic'][
            'child_hom'] = copy.deepcopy(aff_matrix)
        inheritance_report['autosomal']['monoallelic'] = {}
        inheritance_report['autosomal']['monoallelic'][
            'child_het'] = copy.deepcopy(aff_matrix)
        inheritance_report['autosomal']['monoallelic'][
            'child_hom'] = copy.deepcopy(aff_matrix)
        inheritance_report['autosomal']['mosaic'] = {}
        inheritance_report['autosomal']['mosaic']['child_het'] = copy.deepcopy(
            aff_matrix)
        inheritance_report['autosomal']['mosaic']['child_hom'] = copy.deepcopy(
            aff_matrix)
        inheritance_report['autosomal']['imprinted'] = {}
        inheritance_report['autosomal']['imprinted'][
            'child_het'] = copy.deepcopy(aff_matrix)
        inheritance_report['autosomal']['imprinted'][
            'child_hom'] = copy.deepcopy(aff_matrix)

        inheritance_report['allosomal']['hemizygous'] = {}
        inheritance_report['allosomal']['hemizygous'][
            'child_het'] = copy.deepcopy(aff_matrix)
        inheritance_report['allosomal']['hemizygous'][
            'child_hemi'] = copy.deepcopy(aff_matrix)
        inheritance_report['allosomal']['hemizygous'][
            'child_hom'] = copy.deepcopy(aff_matrix)
        inheritance_report['allosomal']['X-linked_dominant'] = {}
        inheritance_report['allosomal']['X-linked_dominant'][
            'child_het'] = copy.deepcopy(aff_matrix)
        inheritance_report['allosomal']['X-linked_dominant'][
            'child_hemi'] = copy.deepcopy(aff_matrix)
        inheritance_report['allosomal']['X-linked_dominant'][
            'child_hom'] = copy.deepcopy(aff_matrix)
        inheritance_report['allosomal']['X-linked_over_dominance'] = {}
        inheritance_report['allosomal']['X-linked_over_dominance'][
            'child_het'] = copy.deepcopy(aff_matrix)
        inheritance_report['allosomal']['X-linked_over_dominance'][
            'child_hemi'] = copy.deepcopy(aff_matrix)
        inheritance_report['allosomal']['X-linked_over_dominance'][
            'child_hom'] = copy.deepcopy(aff_matrix)
        inheritance_report['allosomal']['monoallelic_Y_hemizygous'] = {}
        inheritance_report['allosomal']['monoallelic_Y_hemizygous'][
            'child_hemi'] = copy.deepcopy(aff_matrix)

        return inheritance_report

    def populate_inheritance_report(self, chromtype, mode, child_gt, mum_gt,
                                    dad_gt, mum_aff, dad_aff):
        """
        Populate the inheritance report - both parents present
        """
        if chromtype == 'autosomal':
            if child_gt == '0/1':
                child_geno = 'child_het'
            elif child_gt == '1/1':
                child_geno = 'child_hom'
            else:
                print("Error - unrecognised child GT")
                exit(1)
        elif chromtype == 'allosomal':
            if child_gt == 'hemizygous':
                child_geno = 'child_hemi'
            elif child_gt == 'homozygous':
                child_geno = 'child_hom'
            elif child_gt == 'heterozygous':
                child_geno = 'child_het'
            else:
                print("Error - unrecognised child GT")
                exit(1)
        else:
            print("Error - unrecognised chromosome type")
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
        self.inheritance_report[chromtype][mode][child_geno][dad_aff_state][
            mum_aff_state][parent_gts] += 1
