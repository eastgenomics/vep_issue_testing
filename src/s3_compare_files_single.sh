#!/bin/bash

# # CASE 1 SAMPLE 1 - ARCHIVED!!!
# case_id="MB7226"
# vcf_name="MB7226-226908-1-DNA-egg6_withLowSupportHotspots_annotated.vcf.split.vcf.gz"
# vcf_id="project-GVj5ZF84f2F1K3668FFg0gzq:file-GVpkb7j4g8X7bBjVbB1kb2yx"
# tsv_name="MB7226_CombinedVariantOutput.tsv"
# tsv_id="project-GVj5ZF84f2F1K3668FFg0gzq:file-GVjQX084BGgb2F2J7bp1QVy3"

# # CASE 1 SAMPLE 2
# case_id="123862478-23179S0030-23TSOD21-8471"
# vcf_name="123862478-23179S0030-23TSOD21-8471_withLowSupportHotspots_annotated.vcf.split.vcf.gz"
# vcf_id="project-GY4GKy04V549j9jfQ2v95bp0:file-GY51Jv847b73q2XPy9y6jp5Q"
# tsv_name="123862478-23179S0030-23TSOD21-8471_CombinedVariantOutput.tsv"
# tsv_id="project-GY4GKy04V549j9jfQ2v95bp0:file-GY4g2VQ44FYy0y0pyPzkvq36"

# # CASE 2
# case_id="125125107-23248S0007-23TSOD45-8471"
# vcf_name="125125107-23248S0007-23TSOD45-8471_withLowSupportHotspots_annotated.vcf.split.vcf.gz"
# vcf_id="project-GZ220Qj4zKJ0gZ0JXkQk839z:file-GZ3g1zj47BJq6gqb10Fj1V4Q"
# tsv_name="125125107-23248S0007-23TSOD45-8471_CombinedVariantOutput.tsv"
# tsv_id="project-GZ220Qj4zKJ0gZ0JXkQk839z:file-GZ2Fp6Q456Fj8fQVF7QQg7vX"

# case which ends up with blank p.
case_id="123399008-23156S0026-23TSOD4-8471"
vcf_name="123399008-23156S0026-23TSOD4-8471_withLowSupportHotspots_annotated.vcf.split.vcf.gz"
vcf_id="project-GX8byK047b3PF6KGBGq2b9kx:file-GX9fyg84XQQ9XVQjJzY504V7"
tsv_name="123399008-23156S0026-23TSOD4-8471_CombinedVariantOutput.tsv"
tsv_id="project-GX8byK047b3PF6KGBGq2b9kx:file-GX9Yy1j4X85P6VP425y4P1V0"


# file path shortcuts
tsv_path="${case_id}/${tsv_name}"
vcf_gz_path="${case_id}/${vcf_name}"
vcf_path="${vcf_gz_path%.gz}"
tsv_variants="${case_id}/${case_id}_tsv_variants_v2"

printf "case_id\tgene\tchrom\tpos\tref\talt\ttsv_trscpt\ttsv_pdot\ttsv_vaf\ttsv_cons\tvcf_trscpt_c\tvcf_trscpt_p\tvcf_pdot\tvcf_gnomad_g\tvcf_gnomad_e\tvcf_clnsig\n" \
> "cvo_vs_vcf_output_${case_id}_v2"

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
# i.e. lines where f2 starts with 'chr', f3 is an integer, and f8 isn't blank
# take gene, chrom, pos, ref, alt, vaf, p-dot, consequence
awk -F"\t" '$2~/^chr/ && $3~/^[0-9]+$/ && $8~"p." \
{print $1 "\t" $2 "\t" $3 "\t"  $4 "\t"  $5 "\t" $6 "\t" $8 "\t" $10}' \
"$tsv_path" > "$tsv_variants"

# for each tsv variant, get matching variants from vcf
while read -r gene chrom pos ref alt tsv_vaf trscpt_pdot tsv_cons; do

    # tsv uses e.g. 'chr1' but vcf uses '1'
    chrom="${chrom##chr}"

    # get tsv variant's transcript and p., remove parentheses
    tsv_trscpt=$(echo "$trscpt_pdot" | cut -d ":" -f 1)

    pdot_str=$(echo "$trscpt_pdot" | cut -d ":" -f 2 | grep -Po 'p\.\(.*\)$')
    no_br="${pdot_str//\)/}"
    tsv_pdot="${no_br//\(/}"

    if [[ -n "${tsv_pdot// /}" ]]; then

        # get INFO from vcf lines with the same chrom, pos, ref & alt
        vcf_vars=$(awk -F"\t" \
        -v chrom="$chrom" -v pos="$pos" -v ref="$ref" -v alt="$alt" \
        '$1==chrom && $2==pos && $4==ref && $5==alt {print $8}' "$vcf_path")

        # check whether p. notation is different to the tsv variant
        for var in $vcf_vars; do
            if [[ -n "${var// /}" ]]; then

                # issue with unicode '%3D' not converting to '='
                vcf_var="${var//\%3D/=}"
                vcf_pdot=$(grep -Po '(?<=:)p.(.*?)(?=;)' <<< "$vcf_var")

                # if vcf p. is non-blank and different to tsv p.
                if [[ -n "${vcf_pdot// /}" ]] && \
                [[ "$tsv_pdot" != "$vcf_pdot" ]]; then

                    # get vcf transcripts, gnomad, and clinvar info
                    vcf_trscpt_p=$(grep -Po '(?<=;CSQ_HGVSp=)(.*?)(?=:p.)' <<< "$vcf_var")
                    vcf_trscpt_c=$(grep -Po '(?<=;CSQ_HGVSc=)(.*?)(?=:)' <<< "$vcf_var")

                    vcf_gnomad_g=$(grep -Po '(?<=;CSQ_gnomADg_AF=)(.*?)(?=;)' <<< "$vcf_var")
                    vcf_gnomad_e=$(grep -Po '(?<=;CSQ_gnomADe_AF=)(.*?)(?=;)' <<< "$vcf_var")

                    vcf_clnsig=$(grep -Po '(?<=;CSQ_ClinVar_CLNSIG=)(.*?)(?=;)' <<< "$vcf_var")

                    # if tsv and vcf variants affect the same transcript
                    if [[ "$tsv_trscpt" == "$vcf_trscpt_p" ]] \
                    || [[ "$tsv_trscpt" == "$vcf_trscpt_c" ]]; then

                        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                        "$case_id" "$gene" "$chrom" "$pos" "$ref" "$alt" \
                        "$tsv_trscpt" "$tsv_pdot" "$tsv_vaf" "$tsv_cons" \
                        "$vcf_trscpt_c" "$vcf_trscpt_p" "$vcf_pdot" \
                        "$vcf_gnomad_g" "$vcf_gnomad_e" "$vcf_clnsig" \
                        >> "cvo_vs_vcf_output_${case_id}_v2"
                    fi
                fi
            fi
        done
    fi
done < "$tsv_variants"
