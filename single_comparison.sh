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

echo "Comparing VCF and TSV files for ${case_id}"

# make case dir
if [[ ! -d "$case_id" ]]; then
    mkdir "$case_id"
fi

# download tsv
if [[ ! -e "$tsv_path" ]]; then
    dx download "$tsv_id" -o "$case_id" -f --no-progress
fi

# download and unzip vcf
if [[ ! -e "$vcf_gz_path" ]] && [[ ! -e "$vcf_path" ]]; then
    dx download "$vcf_id" -o "$case_id" -f --no-progress
    gunzip "$vcf_gz_path"
elif [[ -e "$vcf_gz_path" ]]; then
    gunzip "$vcf_gz_path"
fi

# get variant lines from tsv
# i.e. lines where f2 starts with 'chr', f3 is an integer, and f8 contains p.
# take gene, chrom, pos, ref, alt, vaf, p-dot, consequence

tsv_variants=$(awk -F"\t" \
'$2~/^chr/ && $3~/^[0-9]+$/ && $8~"p." \
{print $1 "\t" $2 "\t" $3 "\t"  $4 "\t"  $5 "\t" $6 "\t" $8 "\t" $10}' \
"$tsv_path")

# for each tsv variant, get matching variants from vcf
oldIFS=$IFS
IFS=$'\n'

for variant in $tsv_variants; do

    # read in and assign individual values
    IFS=$'\t' read -r gene chrom pos ref alt tsv_vaf trs_pdot tsv_cons <<< "$variant"

    # tsv uses e.g. 'chr1' but vcf uses '1'
    chrom="${chrom##chr}"

    # get tsv variant's transcript and p., remove parentheses
    tsv_trs=$(echo "$trs_pdot" | cut -d ":" -f 1)

    pdot_str=$(echo "$trs_pdot" | cut -d ":" -f 2 | grep -Po 'p\.\(.*\)$')
    no_br="${pdot_str//\)/}"
    tsv_pdot="${no_br//\(/}"

    # if the tsv's p. value isn't blank,
    if [[ -n "${tsv_pdot// /}" ]]; then

        # get INFO from any vcf lines with the same chrom, pos, ref & alt
        vcf_vars=$(awk -F"\t" \
        -v chrom="$chrom" -v pos="$pos" -v ref="$ref" -v alt="$alt" \
        '$1==chrom && $2==pos && $4==ref && $5==alt {print $8}' "$vcf_path")

        # check whether vcf p. is different to in the tsv
        for var in $vcf_vars; do
            if [[ -n "${var// /}" ]]; then

                # issue with unicode '%3D' not converting to '='
                vcf_var="${var//\%3D/=}"
                vcf_pdot=$(grep -Po '(?<=:)p.(.*?)(?=;)' <<< "$vcf_var")

                # if vcf p. is non-blank and also different to tsv p.
                if [[ -n "${vcf_pdot// /}" ]] && \
                [[ "$tsv_pdot" != "$vcf_pdot" ]]; then

                    # get vcf transcripts, gnomad, and clinvar info
                    vcf_trs_p=$(grep -Po '(?<=;CSQ_HGVSp=)(.*?)(?=:p.)' <<< "$vcf_var")
                    vcf_trs_c=$(grep -Po '(?<=;CSQ_HGVSc=)(.*?)(?=:)' <<< "$vcf_var")

                    vcf_gnomad_g=$(grep -Po '(?<=;CSQ_gnomADg_AF=)(.*?)(?=;)' <<< "$vcf_var")
                    vcf_gnomad_e=$(grep -Po '(?<=;CSQ_gnomADe_AF=)(.*?)(?=;)' <<< "$vcf_var")

                    vcf_clnsig=$(grep -Po '(?<=;CSQ_ClinVar_CLNSIG=)(.*?)(?=;)' <<< "$vcf_var")

                    # if tsv and vcf variants affect the same transcript
                    if [[ "$tsv_trs" == "$vcf_trs_p" ]] \
                    || [[ "$tsv_trs" == "$vcf_trs_c" ]]; then

                        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                        "$case_id" "$gene" "$chrom" "$pos" "$ref" "$alt" \
                        "$tsv_trs" "$tsv_pdot" "$tsv_vaf" "$tsv_cons" \
                        "$vcf_trs_c" "$vcf_trs_p" "$vcf_pdot" \
                        "$vcf_gnomad_g" "$vcf_gnomad_e" "$vcf_clnsig" \
                        >> cvo_vs_vcf_output.tsv
                    fi
                fi
            fi
        done
    fi
done

IFS=$oldIFS
