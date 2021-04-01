"""
Copyright (c) 2021 Genome Research Ltd.
Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import logging
from datetime import datetime

from utils.parse_args import get_options
from file_loading.ped_files import create_ped, openped
from filtering.filter import filter_trio
from output.print_results import create_output

def main():
    """
    run the clinical filtering analyses
    """
    args = get_options()
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    # logfile=""
    if args.ped is None:
        logfile = args.outdir + "/clinical_filter." + timestamp + ".log"
    else:
        logfile = args.ped + "." + timestamp + ".log"
    logging.basicConfig(filename=logfile, level=logging.DEBUG)

    if args.ped is None:
        # create ped file from family members on command line
        args.ped = args.outdir + "/ped." + timestamp + ".ped"
        create_ped(args.ped, args.child, args.mother, args.father, args.sex,
                   args.mum_aff, args.dad_aff)

    families = openped(args.ped, args.proband_list)
    # exit(0)
    variants_per_family = {}
    inheritance_reports_per_family = {}
    for family in families.keys():
        filtered_variants, inheritance_report = filter_trio(families[family], args.known_genes, args.known_regions,
                    args.trusted_variants, args.outdir)
        variants_per_family[family] = filtered_variants
        inheritance_reports_per_family[family] = inheritance_report

    # print(variants_per_family)
    create_output(families, variants_per_family, inheritance_reports_per_family, args.outdir)
    exit(0)
    # vcf_file = '/Volumes/team29/re3/new_clinical_filtering/test_data/vcfs/DDDP100001.gatk.vep.revel.2020-10-29.vcf.gz'
    # readvcf(vcf_file)


if __name__ == "__main__":
    main()
