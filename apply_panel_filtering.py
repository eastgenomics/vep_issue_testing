#!/usr.bin/env python

"""
For each mismatch instance, identifies whether it affects a gene
which is actually in the panel applied to that case.

Requires two other files - a list of the genes in each panel, and
a list of the panels used for each case (from Epic)
"""


import sys

output = []
panel_gene_map = {}
case_panel_map = {}

# read in args (names of input and output files)
input_file = sys.argv[1]
output_file = sys.argv[2]

# read in the panel genes file
# each line is a panel name followed by list of genes

with open('gene_panels_map.txt', 'r') as reader:
    panel_lines = [line.strip() for line in reader.readlines()]

for line in panel_lines:
    split = [val.strip() for val in line.split("\t")]
    panel = split[0]
    panel_gene_map[panel] = split[1:]

# read in the case panels files
# cases may have 1 or 2 panels, so needs to be a list

with open('case_phenotypes.tsv', 'r') as reader:
    case_lines = [line.strip() for line in reader.readlines()]

for line in case_lines:
    split = [val.strip() for val in line.split("\t")]
    case = split[0]
    case_panel_map[case] = split[1:]

# read in all mismatch instances
with open(input_file, 'r') as reader:
    mismatch_lines = [line.strip() for line in reader.readlines()]

# for each mismatch, check whether affected gene is in case's panel(s)
for line in mismatch_lines[1:]:

    # extract relevant variant info
    split = [val.strip() for val in line.split("\t")]

    case, affected_gene, chrom, pos, ref, alt = split[:6]
    tsv_pdot, tsv_vaf = split[7:9]
    vcf_pdot = split[12]
    vcf_gnomADg, vcf_gnomADe = split[14:16]
    variant = f"{chrom}:{pos}:{ref}:{alt}"

    has_no_phenotype = False
    panels_not_listed = False
    no_genes_identified = False
    mismatch_problem = False

    # check case has panel(s) listed
    if case not in case_panel_map.keys():
        has_no_phenotype = True
        case_panel_map[case] = []

    case_panels = case_panel_map[case]

    # identify the genes in those panels
    case_genes = []

    for panel in case_panels:

        if panel not in panel_gene_map.keys():
            panels_not_listed = True
            panel_gene_map[panel] = []

        case_genes += panel_gene_map[panel]

    sorted_genes = sorted(set(case_genes))

    if not case_genes:
        no_genes_identified = True

    # check whether the variant affects a gene in the case's panels
    if affected_gene in case_genes:
        output_dict['mismatch_problem'] = True

    # create output dict
    output_dict = {
        'case': case,
        'panels': case_panels,
        'genes': sorted_genes,
        'variant': variant,
        'affected_gene': affected_gene,
        'tsv_pdot': tsv_pdot,
        'vcf_pdot': vcf_pdot,
        'tsv_vaf': tsv_vaf,
        'vcf_gnomADg': vcf_gnomADg,
        'vcf_gnomADe': vcf_gnomADe,
        'has_no_phenotype': has_no_phenotype,
        'panels_not_listed': panels_not_listed,
        'no_genes_identified': no_genes_identified,
        'mismatch_problem': mismatch_problem}

    output.append(output_dict)

# initialise output file with header line
with open(output_file, 'w') as writer:
    writer.write(
        "case\tpanels\tgenes\tvariant\taffected_gene\t"
        "tsv_pdot\tvcf_pdot\ttsv_vaf\tvcf_gnomADg\tvcf_gnomADe\t"
        "has_no_phenotype\tpanels_not_listed\tno_genes_identified\t"
        "mismatch_problem\n")

# add info for each variant to output
with open(output_file, 'a') as writer:
    for var in output:
        writer.write(
            f"{var['case']}\t{var['panels']}\t{var['genes']}\t"
            f"{var['variant']}\t{var['affected_gene']}\t"
            f"{var['tsv_pdot']}\t{var['vcf_pdot']}\t{var['tsv_vaf']}\t"
            f"{var['vcf_gnomADg']}\t{var['vcf_gnomADe']}\t"
            f"{var['has_no_phenotype']}\t{var['panels_not_listed']}\t"
            f"{var['no_genes_identified']}\t{var['mismatch_problem']}\n")
