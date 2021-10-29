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


def load_genes(genes_file):
    """
    load genes from DDG2P
    """
    genes = {}
    with open(genes_file, 'r') as g:
        lines = g.readlines()
        for l in lines:
            if not l.startswith('chr'):
                linedata = (l).split("\t")
                chr = linedata[0]
                start = linedata[1]
                end = linedata[2]
                symbol = linedata[3]
                hgnc_id = linedata[4]
                status = linedata[5]
                mode = linedata[6]
                mechanism = linedata[7]
                if hgnc_id in genes.keys():
                    genes[hgnc_id]['status'].add(status)
                    genes[hgnc_id]['mode'].add(mode)
                    genes[hgnc_id]['mechanism'].add(mechanism)
                else:
                    genes[hgnc_id] = {}
                    genes[hgnc_id]['chr'] = chr
                    genes[hgnc_id]['start'] = start
                    genes[hgnc_id]['end'] = end
                    genes[hgnc_id]['symbol'] = symbol
                    genes[hgnc_id]['status'] = set({status})
                    genes[hgnc_id]['mode'] = set({mode})
                    genes[hgnc_id]['mechanism'] = set({mechanism})

    return genes


def load_regions():
    """
    load regions of interest
    """
    pass


def load_trusted_variants():
    """
    load trusted variants list - eg clinvar
    """
    pass
