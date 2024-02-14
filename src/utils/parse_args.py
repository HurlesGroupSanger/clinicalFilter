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

import argparse
import os


def get_options():
    """
    Get the options from the command line
    """
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--ped",
        help="Path to ped file containing cohort details for " "multiple trios.",
    )
    group.add_argument("--child", help="Path to child's VCF file.")
    parser.add_argument("--sex", help="The child's chromosomal sex (eg XY, XX, XXY).")
    parser.add_argument("--mother", help="Path to mother's VCF file.")
    parser.add_argument("--father", help="Path to father's VCF file.")
    parser.add_argument("--mum-aff", help="Mother's affected status (1=unaffected, or " "2=affected).")
    parser.add_argument("--dad-aff", help="Father's affected status (1=unaffected, or " "2=affected).")
    parser.add_argument("--proband-list", help="List of probands to be analysed.")
    parser.add_argument("--known-genes", help="Path to file of known disease causative genes.")

    parser.add_argument("--known-regions", help="Path to file of known disease causative regions.")

    parser.add_argument(
        "--trusted-variants",
        help="Path to file of known disease causative " "variants.",
    )

    parser.add_argument("--outdir", help="Output directory.")

    args = parser.parse_args()

    if args.child is not None:
        if args.father is not None and args.dad_aff is None:
            parser.error("--dad-aff must also be used if --father is used")
        if args.mother is not None and args.mum_aff is None:
            parser.error("--mum-aff must also be used if --mother is used")
        if args.sex is None:
            parser.error("--sex must also be used if --child is used")

    if args.outdir is None:
        args.outdir = os.getcwd()

    return args
