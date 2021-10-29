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

import logging


class PostInheritanceFiltering(object):
    """
    Post-inheritance filters
    """

    def __init__(self, candidate_variants, family):
        self.candidate_variants = candidate_variants
        self.family = family

    def postinheritance_filter(self):
        """
        Post-inheritance filtering - MAF and REVEL
        """
        self.maf_filter()
        self.allele_count_filter()
        return self.candidate_variants

    def maf_filter(self):
        """
        Filter non-Biallelic vairants with more stringent MAF thresholds
        """
        for v in list(self.candidate_variants['single_variants'].keys()):
            if self.candidate_variants['single_variants'][v][
                'mode'] != 'Biallelic':
                ddd_af = self.candidate_variants['single_variants'][v][
                    'variant'].ddd_af
                max_af = self.candidate_variants['single_variants'][v][
                    'variant'].max_af
                if ddd_af == '.':
                    ddd_af = '0'
                if max_af == '.':
                    max_af = '0'
                maximum_af = max(float(ddd_af), float(max_af))

                if self.family.has_both_parents() and maximum_af >= 0.0005:
                    del self.candidate_variants['single_variants'][v]
                    logging.info(
                        v + " failed post-inhertance MAF filter for family "
                            "with parents, max AF = " + str(
                            maximum_af))
                elif not self.family.has_both_parents() and maximum_af >= 0.0001:
                    del self.candidate_variants['single_variants'][v]
                    logging.info(
                        v + " failed post-inhertance MAF filter for family "
                            "without parents, max AF = " + str(
                            maximum_af))

    def allele_count_filter(self):
        """
        Filter on AC_het and AC_hemi
        """
        for v in list(self.candidate_variants['single_variants'].keys()):
            if self.candidate_variants['single_variants'][v][
                'mode'] == 'Monoallelic':
                if int(self.candidate_variants['single_variants'][v][
                           'variant'].AC_het) > 4:
                    del self.candidate_variants['single_variants'][v]
                    logging.info(
                        v + " failed post-inhertance AC_het filter for "
                            "monoallelic genes " +
                        self.candidate_variants['single_variants'][v][
                            'variant'].AC_het)
            if self.candidate_variants['single_variants'][v][
                'mode'] == 'Hemizygous' and \
                    self.candidate_variants['single_variants'][v][
                        'sex'] == "XY":
                if int(self.candidate_variants['single_variants'][v][
                           'variant'].AC_hemi) > 0:
                    del self.candidate_variants['single_variants'][v]
                    logging.info(
                        v + " failed post-inhertance AC_hemi filter for "
                            "monoallelic genes " +
                        self.candidate_variants['single_variants'][v][
                            'variant'].AC_hemi)
            if self.candidate_variants['single_variants'][v][
                'mode'] == 'X-linked dominant':
                AC_total = int(self.candidate_variants['single_variants'][v][
                                   'variant'].AC_het) + int(
                    self.candidate_variants['single_variants'][v][
                        'variant'].AC_hemi)
                if AC_total > 4:
                    del self.candidate_variants['single_variants'][v]
                    logging.info(
                        v + " failed post-inhertance AC_hemi filter for "
                            "X linked dominant genes " +
                        self.candidate_variants['single_variants'][v][
                            'variant'].AC_het + " + " +
                        self.candidate_variants['single_variants'][v][
                            'variant'].AC_hemi)
