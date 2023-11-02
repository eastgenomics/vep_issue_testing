#!/bin/bash

# read in args (passed by get_file_info.sh)
case_id="$1"
vcf_id="$2"
vcf_name="$3"
tsv_id="$4"
tsv_name="$5"

tsv_path="${case_id}/${tsv_name}"
vcf_gz_path="${case_id}/${vcf_name}"
vcf_path="${vcf_gz_path%.gz}"

# make case dir
if [[ ! -d "$case_id" ]]; then
    mkdir "$case_id"
fi

# download tsv
if [[ ! -e "$tsv_path" ]]; then
    dx download "$tsv_id" -o "$case_id" -f --no-progress
fi

# download vcf
if [[ ! -e "$vcf_gz_path" ]] && [[ ! -e "$vcf_path" ]]; then
    dx download "$vcf_id" -o "$case_id" -f --no-progress
fi

# keep tsv lines where f2 starts 'chr', f3 is integer, and f8 contains p.
# but not if VAF < 0.05 (one of the workbooks filters)
# take gene, chrom, pos, ref, alt, vaf, p-dot, type

tsv_variants=$(awk -F"\t" \
'$2~/^chr/ && $3~/^[0-9]+$/ && $8~"p." && $6>=0.05 \
{print $1 "\t" $2 "\t" $3 "\t"  $4 "\t"  $5 "\t" $6 "\t" $8 "\t" $10}' \
"$tsv_path")

# for each tsv variant, get matching variants from vcf
oldIFS=$IFS
IFS=$'\n'

for variant in $tsv_variants; do

    # read in and assign individual values
    IFS=$'\t' read -r gene chrom pos ref alt tsv_vaf trs_pdot tsv_type <<< "$variant"

    # tsv uses e.g. 'chr1' but vcf uses '1'
    chrom="${chrom##chr}"

    # get tsv variant's transcript and p., remove parentheses
    tsv_trs=$(echo "$trs_pdot" | cut -d ":" -f 1)

    pdot_str=$(echo "$trs_pdot" | cut -d ":" -f 2 | grep -Po 'p\.\(.*\)$')
    no_br="${pdot_str//\)/}"
    tsv_pdot="${no_br//\(/}"

    # if tsv p. isn't blank,
    if [[ -n "${tsv_pdot// /}" ]]; then

        # make sure vcf is bgzipped (not gzipped), then index
        if [[ -e "${vcf_path}.gz" ]]; then
            gunzip "$vcf_path"
        fi
        if [[ -e "$vcf_path" ]]; then
            bgzip "$vcf_path" -f
        fi
        bcftools index "${vcf_path}.gz" -f

        # retrieve INFO fields where chrom, pos, ref, alt and transcript match
        vcf_fields="$(bcftools query \
        -f '%INFO/DP\t%CSQ_SYMBOL\t%CSQ_HGVSc\t%CSQ_HGVSp\t%CSQ_gnomADg_AF\t%CSQ_gnomADe_AF\t%CSQ_Consequence\n' \
        -r "${chrom}:${pos}" \
        -i "REF='${ref}' && ALT='${alt}' && (INFO/CSQ_HGVSc~'${tsv_trs}' || INFO/CSQ_HGVSp~'${tsv_trs}')" \
        "${vcf_path}.gz")"

        # assign INFO fields to variables
        IFS=$'\t' read -r vcf_dp vcf_symbol vcf_HGVSc vcf_HGVSp vcf_gnomADg vcf_gnomADe vcf_consq <<< "$vcf_fields"

        # parse out the vcf transcripts and p.
        vcf_trs_c=$(echo "$vcf_HGVSc" | grep -Po '^(.*?)(?=:c\.)')
        vcf_trs_p=$(echo "$vcf_HGVSp" | grep -Po '^(.*?)(?=:p\.)')
        vcf_pdot=$(echo "$vcf_HGVSp" | grep -Po '(?<=:)p\.(.*?)$')

        # if the vcf pdot doesn't match,
        if [[ "$tsv_pdot" != "$vcf_pdot" ]]; then

            # list consequences to exclude variants with
            ignore_cons="intron_variant&non_coding_transcript_variant non_coding_transcript_exon_variant 3_prime_UTR_variant 5_prime_UTR_variant downstream_gene_variant intron_variant splice_region_variant&intron_variant splice_region_variant&synonymous_variant synonymous_variant"

            # compare variant to filters
            keep_var=true

            if [[ "$vcf_dp" -eq 0 ]] \
            || [[ $(bc -l <<< "$vcf_gnomADe > 0.01") -eq 1 ]] \
            || [[ $(bc -l <<< "$vcf_gnomADg > 0.01") -eq 1 ]] \
            || [[ "$ignore_cons" == *"$vcf_consq"* ]]; then

                keep_var=false
            fi

            if [[ "$vcf_consq" == "upstream_gene_variant" ]] \
            && [[ "$vcf_symbol" != "TERT" ]]; then

                keep_var=false
            fi

            # if variant passes filters, add to output
            if [[ "$keep_var" == true ]]; then

                printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                "$case_id" "$gene" "$chrom" "$pos" "$ref" "$alt" \
                "$tsv_trs" "$tsv_pdot" "$tsv_vaf" "$tsv_type" \
                "$vcf_trs_c" "$vcf_trs_p" "$vcf_pdot" "$vcf_dp" \
                "$vcf_gnomADg" "$vcf_gnomADe" "$vcf_consq" \
                >> cvo_vs_vcf_output.tsv
            fi
        fi
    fi
done

IFS=$oldIFS
