"""
copyright
"""
import subprocess
import logging
import os

from variants.snv import SNV
from variants.cnv import CNV


def load_variants(family, outdir, regions=None):
    '''get variants in child and parents'''
    proband_vcf = family.proband.get_vcf_path()

    regionfile = outdir + "/reg.tmp"
    with open(regionfile, 'w') as rf:
        for r in regions:
            rf.write(r + "\n")

    sortedregs = outdir + "/reg.tmp_sorted"
    sortcmd = "sort -k1,1V -k2,2n -k3,3n " + regionfile + " > " + sortedregs
    os.system(sortcmd)

    child_vars = readvcf(proband_vcf, sortedregs, family.proband.get_sex())

    mum_vars = {}
    dad_vars = {}
    if not family.has_no_parents():
        # get a region string from the child variants as we don't need parental
        # variants which are not in the child

        childregionfile = outdir + '/childreg.tmp'
        with open(childregionfile, 'w') as crf:
            for varid in child_vars.keys():
                idsplit = varid.split("_")
                reg = idsplit[0] + "\t" + idsplit[1]
                crf.write(reg + "\n")

        childsortedregs = outdir + "/childreg.tmp_sorted"
        childsortcmd = "sort -k1,1V -k2,2n " + childregionfile + " > " + childsortedregs

        os.system(childsortcmd)

        if family.has_mum():
            mum_vcf = family.mum.get_vcf_path()
            mum_vars = readvcf(mum_vcf, childsortedregs, 'F')

        if family.has_dad():
            dad_vcf = family.dad.get_vcf_path()
            dad_vars = readvcf(dad_vcf, childsortedregs, 'M')

        # remove child regions files used for bcftools queries
        os.system("rm " + childregionfile)
        os.system("rm " + childsortedregs)

    # remove regions files used for bcftools queries
    os.system("rm " + regionfile)
    os.system("rm " + sortedregs)

    variants = {}
    variants['child'] = child_vars
    variants['mum'] = mum_vars
    variants['dad'] = dad_vars
    return variants


def readvcf(filename, regions, sex):
    """
    read vcf files and return a dict of variant objects
    """
    vars = {}

    # get list of info fields in the vcf
    fieldlistcmd = "bcftools view -h " + filename + " z | grep ^##INFO | sed 's/^.*ID=// ; s/,.*//'"
    fieldlist = runcommand(fieldlistcmd)
    fields = fieldlist.split("\n")

    # infofields_wanted = ['Consequence', 'Gene', 'SYMBOL', 'Feature', 'CANONICAL',
    #               'MANE', 'HGNC_ID', 'MAX_AF', 'MAX_AF_POPS',
    #               'DDD_AF', 'REVEL', 'PolyPhen', 'Protein_position', 'HGVSc', 'HGVSp']
    # formatfields = ['GT', 'GQ', 'PID', 'AD']

    infofields_wanted = ['Consequence', 'Gene', 'SYMBOL', 'Feature',
                         'CANONICAL',
                         'MANE', 'HGNC_ID', 'MAX_AF', 'MAX_AF_POPS',
                         'DDD_AF', 'DDD_father_AF', 'REVEL', 'PolyPhen',
                         'Protein_position',
                         'HGVSc', 'HGVSp', 'DNM_TYPE', 'END', 'SVTYPE', 'SVLEN',
                         'CNVFILTER', 'HGNC_ID_ALL']
    formatfields = ['GT', 'GQ', 'PID', 'AD', 'CIFER_INHERITANCE', 'CN']

    # create infostring containing only the fields present
    info_query = []
    for inf in infofields_wanted:
        info_query.append("%INFO/" + inf)
        # if inf in fields:
        #     info_query.append("%INFO/" + inf)
        # else:
        #     info_query.append("")

    infostring = ("\t").join(info_query)
    formatstring = ("\t%").join(formatfields)

    if regions is None:
        bcfcmdroot = "bcftools norm -m - " + filename + \
                     " | bcftools view -e 'INFO/MAX_AF>0.005 | FORMAT/GT[0]!=" + \
                     '"alt"' + "'  | bcftools query -u -f '%CHROM\t%POS\t%REF\t%ALT{0}\t"
    else:
        bcfcmdroot = "bcftools norm -m - -R " + regions + " " + filename + \
                     " | bcftools view -e 'INFO/MAX_AF>0.005 | FORMAT/GT[0]!=" + \
                     '"alt"' + "'  | bcftools query -u -f '%CHROM\t%POS\t%REF\t%ALT{0}\t"

    bcfcmd = bcfcmdroot + infostring + "[\t%" + formatstring + "]\n'"
    output = runcommand(bcfcmd)

    if not output == "Error in command":
        outputlines = output.split('\n')
    else:
        logging.error("Variants not loaded from " + filename)

    for ol in outputlines:
        oldata = ol.split("\t")
        alt = oldata[3]
        if alt == '*':  # get rid of any where alt allele is *
            continue
        # populate hash with variant data
        varid = ("_").join([oldata[0], oldata[1], oldata[2], alt])
        vdata = {}
        vdata['chrom'] = oldata[0]
        vdata['pos'] = oldata[1]
        vdata['ref'] = oldata[2]
        vdata['alt'] = alt
        vdata['consequence'] = oldata[4]
        vdata['ensg'] = oldata[5]
        vdata['symbol'] = oldata[6]
        vdata['feature'] = oldata[7]
        vdata['canonical'] = oldata[8]
        vdata['mane'] = oldata[9]
        vdata['hgnc_id'] = oldata[10]
        vdata['max_af'] = oldata[11]
        vdata['max_af_pops'] = oldata[12]
        vdata['ddd_af'] = oldata[13]
        vdata['ddd_father_af'] = oldata[14]
        vdata['revel'] = oldata[15]
        vdata['polyphen'] = oldata[16]
        vdata['protein_position'] = oldata[17]
        vdata['hgvsc'] = oldata[18]
        vdata['hgvsp'] = oldata[19]
        vdata['sex'] = sex
        vdata['dnm'] = oldata[20]
        vdata['cnv_end'] = oldata[21]
        vdata['cnv_type'] = oldata[22]
        vdata['cnv_length'] = oldata[23]
        vdata['cnv_filter'] = oldata[24]
        vdata['hgnc_id_all'] = oldata[25]
        vdata['gt'] = oldata[26]
        vdata['gq'] = oldata[27]
        vdata['pid'] = oldata[28]
        vdata['ad'] = oldata[29]
        vdata['cnv_inh'] = oldata[30]
        vdata['copy_num'] = oldata[31]

        # if vdata['dnm'] == 'DNM':
        #     if len(vdata['ref']) == len(vdata['alt']):
        #         vdata['denovo_snv'] = True
        #     else:
        #         vdata['denovo_indel'] = True
        # else:
        #     if len(vdata['ref']) == len(vdata['alt']):
        #         vdata['denovo_snv'] = False
        #     else:
        #         vdata['denovo_indel'] = False

        var = SNV
        if alt in ['<DEL>', '<DUP>']:
            var = CNV
        vars[varid] = var(vdata)

    logging.info("Variants loaded from " + filename)

    return vars


def runcommand(cmd):
    try:
        byteOutput = subprocess.check_output(cmd, shell=True)
        return byteOutput.decode('UTF-8').rstrip()
    except subprocess.CalledProcessError as e:
        return "Error in command"
