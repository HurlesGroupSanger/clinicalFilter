"""
copyright
"""

def load_genes(genes_file):
    '''load genes from DDG2P'''
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
    '''load regions of interest'''
    pass

def load_trusted_variants():
    '''load trusted variants list - eg clinvar'''
    pass
