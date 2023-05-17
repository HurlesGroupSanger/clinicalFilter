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


def common_elements(list1, list2):
    return list(set(list1) & set(list2))


def add_single_var_to_candidates(varid, var, hgncid, inh, candidates):
    """
    Add a single variant to candidate variants hash
    """
    if not varid in candidates["single_variants"].keys():
        candidates["single_variants"][varid] = {}
        candidates["single_variants"][varid]["mode"] = set()
        candidates["single_variants"][varid]["mode"].add(inh)
        candidates["single_variants"][varid]["variant"] = var
        candidates["single_variants"][varid]["hgncid"] = hgncid
    else:
        candidates["single_variants"][varid]["mode"].add(inh)


def add_compound_het_to_candidates(varid, var, hgncid, inh, candidates):
    """
    Add a potential compound het to candidate variants hash
    """
    # todo could probably streamline this

    if not hgncid in candidates["compound_hets"].keys():
        candidates["compound_hets"][hgncid] = {}
        candidates["compound_hets"][hgncid][varid] = {}
        candidates["compound_hets"][hgncid][varid]["mode"] = set()
        candidates["compound_hets"][hgncid][varid]["mode"].add(inh)
        candidates["compound_hets"][hgncid][varid]["variant"] = var
        candidates["compound_hets"][hgncid][varid]["hgncid"] = hgncid
    else:
        if not varid in candidates["compound_hets"][hgncid].keys():
            candidates["compound_hets"][hgncid][varid] = {}
            candidates["compound_hets"][hgncid][varid]["mode"] = set()
            candidates["compound_hets"][hgncid][varid]["mode"].add(inh)
            candidates["compound_hets"][hgncid][varid]["variant"] = var
            candidates["compound_hets"][hgncid][varid]["hgncid"] = hgncid
        else:
            candidates["compound_hets"][hgncid][varid]["mode"].add(inh)


def convert_genotype_to_gt(genotype):
    """
    Convert 0/1/2 genotype to GATK genotype
    """
    gt = ""
    if genotype == "0":
        gt = "0/0"
    elif genotype == "1":
        gt = "0/1"
    elif genotype == "2":
        gt = "1/1"
    else:
        print("Invalid genotype - shouldn't get here")
        exit(1)

    return gt
