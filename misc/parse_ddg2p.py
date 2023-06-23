#! /usr/bin/python

#parse information from gene2phenotype and ensembl to create ddg2p file

import os
import re
import csv
import pprint

ddg2pfile = '/lustre/scratch115/projects/ddd/users/ddd/df10kv2_2_clinical_filter/ddg2p/DDG2P_6_8_2020.csv'


ensemblfile = "/software/ddd/resources/ensembl/ensembl_biomart_export_190402.txt"


outfile = '/lustre/scratch115/projects/ddd/users/ddd/df10kv2_2_clinical_filter/ddg2p/DDG2P_20200810_clinical.filtering.txt'

#can have >1 line per gene in the ddg2p file (different syndromes)
#diff syndromes can have diff mechanisms

ddg2p = {}
genes = {}
ensgs = {}

index = 0

wanted_cats = {'both RD and IF':'Both RD and IF', 'confirmed':'Confirmed DD gene', 'probable':'Probable DD gene'}

#get info out of ddg2p and do some basic formatting
with open(ddg2pfile, 'r') as d:
    reader = csv.DictReader(d)
    for row in reader:
        index += 1
        if row['DDD category'] in wanted_cats.keys():
            if row['mutation consequence'] == '':
                continue
            if row['allelic requirement'] == '':
                continue
            ddg2p[index] = {}
            ddg2p[index]['gene'] = row['gene symbol']
            ddg2p[index]['hgnc_id'] = row['hgnc id']
            ddg2p[index]['type'] = wanted_cats[row['DDD category']]#reformat category
            ddg2p[index]['mode'] = row['allelic requirement'].capitalize()
            ddg2p[index]['mechanism'] = row['mutation consequence'].capitalize()
            ddg2p[index]['syndrome'] = row['disease name'] 
            if not row['hgnc id'] in genes.keys():
                genes[row['hgnc id']] = {}
                genes[row['hgnc id']]['symbol'] = row['gene symbol']


chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT']
#get gene co-ordinates from ensembl biomart download
with open(ensemblfile, 'r') as e:
    elines = e.readlines()
    for line in elines:
        ldata = line.split("\t")
        if not ldata[0] == 'Gene stable ID':
            hgnc = ldata[4]
            gene = ldata[3]
            chrom = ldata[5].strip()
            if chrom not in chromosomes:#this gets rid of entries in patches etc so should remove duplicates
                continue
            start = ldata[1]
            end = ldata[2]
            if hgnc in genes.keys():
                if 'chrom' in genes[hgnc].keys():#check for remaining duplicates
                    print ">1 row for hgnc " + hgnc + " " + gene + " in ensembl"
                else:
                    genes[hgnc]['chrom'] = chrom
                    genes[hgnc]['start'] = start
                    genes[hgnc]['end'] = end


header = ['chr', 'start', 'stop', 'gene', 'hgnc_id', 'type', 'mode', 'mech', 'syndrome']

with open(outfile, 'w') as o:
    o.write(("\t").join(header))
    o.write("\n")
    for item in ddg2p.keys():
        hgnc_id = ddg2p[item]['hgnc_id'] 
        linedata = []
        linedata.append(genes[hgnc_id]['chrom'])
        linedata.append(genes[hgnc_id]['start'])
        linedata.append(genes[hgnc_id]['end'])
        linedata.append(ddg2p[item]['gene'])
        linedata.append(hgnc_id)
        linedata.append(ddg2p[item]['type'])
        linedata.append(ddg2p[item]['mode'])
        linedata.append(ddg2p[item]['mechanism'])
        linedata.append(ddg2p[item]['syndrome'])
        o.write(("\t").join(linedata))
        o.write("\n")
