#!/usr/bin/env python


variants = {}

# read in all instances of variant mismatches between tsv and vcf
with open('20231017_cvo_vs_vcf_output_all_v2', 'r') as reader:
    lines = reader.readlines()

for line in lines[1:]:  # ignore header

    case_id, gene, chrom, pos, ref, alt, tsv_trscpt, tsv_pdot, tsv_vaf, \
        tsv_cons, vcf_trscpt_c, vcf_trscpt_p, vcf_pdot, vcf_gnom_g, \
        vcf_gnom_e, vcf_clinvar = line.split('\t')

    # unique variants include position, change, and transcript
    variant = f"{chrom}:{pos}:{ref}:{alt}:{tsv_trscpt}"

    if variant not in variants.keys():

        variants[variant] = {
            'gene': gene.strip(),
            'tsv_pdot': tsv_pdot.strip(),
            'tsv_cons': tsv_cons.strip(),
            'vcf_trscpt_c': vcf_trscpt_c.strip(),
            'vcf_trscpt_p': vcf_trscpt_p.strip(),
            'vcf_pdot': vcf_pdot.strip(),
            'vcf_gnom_g': vcf_gnom_g.strip(),
            'vcf_gnom_e':vcf_gnom_e.strip(),
            'vcf_clinvar': vcf_clinvar.strip(),
            'cases': {case_id.strip(): tsv_vaf.strip()}
        }

    # if variant hasn't been seen, check pdots are the same as already in dict
    else:
        tsv_pdot_seen = variants[variant]['tsv_pdot']
        vcf_pdot_seen = variants[variant]['vcf_pdot']
        vcf_tr_p_seen = variants[variant]['vcf_trscpt_p']
        vcf_tr_c_seen = variants[variant]['vcf_trscpt_c']

        assert tsv_pdot == tsv_pdot_seen; \
            f"multiple different tsv_pdot for {variant}"

        assert vcf_pdot == vcf_pdot_seen; \
            f"multiple different vcf_pdot for {variant}"

        assert vcf_trscpt_p == vcf_tr_p_seen; \
            f"multiple different vcf_trscpt_p for {variant}"

        assert vcf_trscpt_c == vcf_tr_c_seen; \
            f"multiple different vcf_trscpt_c for {variant}"

        variants[variant]['cases'][case_id] = tsv_vaf

unique_var_count = len(variants.keys())
print(f"{unique_var_count} unique variants have mismatching annotation")

with open('unique_mismatch_variants', 'w') as writer:
    writer.write("chrom\tpos\tref\talt\tgene\ttsv_pdot\ttsv_trscpt\ttsv_cons\tvcf_pdot\tvcf_trscpt_p\tvcf_trscpt_c\tvcf_gnomad_g\tvcf_gnomad_e\tvcf_clinvar\tcase_count\tcases\n")

for var, var_info in variants.items():

    chrom, pos, ref, alt, tsv_trscpt = var.split(':')

    case_list = [f"{case_id} ({vaf})" \
        for case_id, vaf in variants[var]['cases'].items()]

    case_count = str(len(case_list))
    case_str = ', '.join(case_list).strip()

    line = '\t'.join(
        [chrom, pos, ref, alt, var_info['gene'],
        var_info['tsv_pdot'], tsv_trscpt, var_info['tsv_cons'],
        var_info['vcf_pdot'], var_info['vcf_trscpt_p'],
        var_info['vcf_trscpt_c'], var_info['vcf_gnom_g'],
        var_info['vcf_gnom_e'], var_info['vcf_clinvar'], case_count, case_str])

    with open('unique_mismatch_variants', 'a') as writer:
        writer.write(f"{line}\n")
