"""
copyright
"""

import argparse
import os

def get_options():
    """gets the options from the command line
    """
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--ped",
                       help="Path to ped file containing cohort details for multiple trios.")
    group.add_argument("--child", help="Path to child's VCF file.")
    parser.add_argument("--gender", help="The child's gender (M or F).")
    parser.add_argument("--mother", help="Path to mother's VCF file.")
    parser.add_argument("--father", help="Path to father's VCF file.")
    parser.add_argument("--mum-aff",
                        help="Mother's affected status (1=unaffected, or 2=affected).")
    parser.add_argument("--dad-aff",
                        help="Father's affected status (1=unaffected, or 2=affected).")
    parser.add_argument("--proband-list", help="List of probands to be analysed.")
    parser.add_argument("--known-genes",
                        help="Path to file of known disease causative genes.")

    parser.add_argument("--known-regions",
                        help="Path to file of known disease causative regions.")

    parser.add_argument("--trusted-variants",
                        help="Path to file of known disease causative variants.")

    parser.add_argument("--outdir", help="Output directory.")

    args = parser.parse_args()

    if args.child is not None:
        if args.father is not None and args.dad_aff is None:
            parser.error("--dad-aff must also be used if --father is used")
        if args.mother is not None and args.mum_aff is None:
            parser.error("--mum-aff must also be used if --mother is used")
        if args.gender is None:
            parser.error("--gender must also be used if --child is used")
        genders = ['M', 'F']
        if not args.gender in genders:
            parser.error("--gender must be M or F")

    if args.outdir is None:
        args.outdir = os.getcwd()

    return args
