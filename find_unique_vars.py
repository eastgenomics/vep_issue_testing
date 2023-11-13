#!/usr/bin/env python

## once a list of all mismatch instances has been generated, processes this to
## list unique mismatches and the specific cases each affects


variants = {}
seen_variants=[]

# read in args (name of input and output files)
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

# read in lines from input file
with open(input_file, 'r') as reader:
    lines = reader.readlines()

for line in lines[1:]:  # ignore header line

    # split line into individual variables
    case_id, gene, chrom, pos, ref, alt, tsv_trs, tsv_pdot, tsv_vaf, \
        tsv_type, vcf_trs_c, vcf_trs_p, vcf_pdot, vcf_dp, vcf_gnomADg, \
        vcf_gnomADe, vcf_consq = line.split('\t')

    # unique variants include position, change, and transcript
    variant = f"{chrom}:{pos}:{ref}:{alt}:{tsv_trs}"

    # if variant hasn't been seen already, create a dict for it
    if variant not in seen_variants:

        seen_variants.append(variant)

        variants[variant] = {
            'gene': gene.strip(),
            'tsv_pdot': tsv_pdot.strip(),
            'tsv_type': tsv_type.strip(),
            'vcf_trs_c': vcf_trs_c.strip(),
            'vcf_trs_p': vcf_trs_p.strip(),
            'vcf_pdot': vcf_pdot.strip(),
            'vcf_gnomADg': vcf_gnomADg.strip(),
            'vcf_gnomADe':vcf_gnomADe.strip(),
            'vcf_consq': vcf_consq.strip(),
            'cases': [case_id.strip()]}

    # if the variant already has a dict,
    else:
        # check tsv/vsf pdots are the same as already seen
        tsv_pdot_seen = variants[variant]['tsv_pdot']
        vcf_pdot_seen = variants[variant]['vcf_pdot']
        vcf_trs_p_seen = variants[variant]['vcf_trs_p']
        vcf_trs_c_seen = variants[variant]['vcf_trs_c']

        assert tsv_pdot == tsv_pdot_seen; \
            f"multiple different tsv_pdot for {variant}"

        assert vcf_pdot == vcf_pdot_seen; \
            f"multiple different vcf_pdot for {variant}"

        assert vcf_trs_p == vcf_trs_p_seen; \
            f"multiple different vcf_trs_p for {variant}"

        assert vcf_trs_c == vcf_trs_c_seen; \
            f"multiple different vcf_trs_c for {variant}"

        # add case id to list in dict
        variants[variant]['cases'].append(case_id.strip())

# write the output file header
with open(output_file, 'w') as writer:
    writer.write("chrom\tpos\tref\talt\tgene\ttsv_pdot\ttsv_trs\ttsv_type\tvcf_pdot\tvcf_trs_p\tvcf_trs_c\tvcf_gnomad_g\tvcf_gnomad_e\tvcf_consq\tcase_count\tcases\n")

# append the info for each unique variant to the output file
for var, var_info in variants.items():

    chrom, pos, ref, alt, tsv_trs = var.split(':')

    # get list of cases affected by that variant
    case_list = variants[var]['cases']
    case_count = str(len(case_list))
    case_str = ', '.join(case_list).strip()

    # create the line to output
    line = '\t'.join(
        [chrom, pos, ref, alt, var_info['gene'],
        var_info['tsv_pdot'], tsv_trs, var_info['tsv_type'],
        var_info['vcf_pdot'], var_info['vcf_trs_p'],
        var_info['vcf_trs_c'], var_info['vcf_gnomADg'],
        var_info['vcf_gnomADe'], var_info['vcf_consq'], case_count, case_str])

    # append it to the file
    with open(output_file, 'a') as writer:
        writer.write(f"{line}\n")
