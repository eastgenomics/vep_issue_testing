#!/bin/bash

# define output files
date=$(date '+%Y%m%d')
output_all_vars="output_all_mismatches_${date}.tsv"
output_unique_vars="output_unique_mismatches_${date}.tsv"
output_panel_filtering="output_panel_filtering_${date}.tsv"

# initialise main output file
printf "case_id\tgene\tchrom\tpos\tref\talt\ttsv_trscpt\ttsv_pdot\t \
tsv_vaf\ttsv_type\tvcf_trscpt_c\tvcf_trscpt_p\tvcf_pdot\tvcf_dp\tvcf_gnomad_g \
\tvcf_gnomad_e\tvcf_consq\n" > "$output_all_vars"

# identify all 002_*_TSO500 projects between 22 May and 20 Oct 2023
projects=$(dx find projects \
--name "002_*_TSO500" \
--created-after 2023-05-20 \
--created-before 2023-10-20 \
--brief)

# find all the relevant files from each project
i=1
for project in $projects; do

    echo "Processing project ${i}: ${project}"

    # get IDs of all VEP-annotated VCF files in that project
    vcf_ids=$(dx find data \
    --name "*withLowSupportHotspots_annotated.vcf.split.vcf.gz" \
    --project "$project" \
    --brief)

    # for each VCF, get the file name
    for vcf_id in $vcf_ids; do

        vcf_name=$(dx describe "$vcf_id" --json | jq -r '.name')
        vcf_prefix=$(echo "$vcf_name" | cut -f 1 -d $'_')

        # extract the case id
        if [[  $vcf_prefix == *"egg6" ]]; then
            case_id=$(echo "$vcf_prefix" | cut -f 1 -d $'-')
        else
            case_id="$vcf_prefix"
        fi

        # ignore QC samples
        if [[ "$case_id" != *"HD753"* ]] \
        && [[ "$case_id" != *"Q"* ]]; then

            # identify the TSO500 output file - same tsv file appears
            # twice in a project, so only keep one

            tsv_id=$(dx find data \
                --name "${case_id}_CombinedVariantOutput.tsv" \
                --project "$project" \
                --brief | cut -d $'\n' -f 1)

            # if there is such a matching tsv file (tsv_id isn't blank),
            if [[ "${tsv_id// /}" != "" ]]; then

                # check whether the vcf or tsv are archived
                vcf_state=$(dx describe "$vcf_id" --json | jq -r '.archivalState')
                tsv_state=$(dx describe "$tsv_id" --json | jq -r '.archivalState')

                # if both files are live,
                if [[ "$vcf_state" == 'live' ]] \
                && [[ "$tsv_state" == 'live' ]]; then

                    # get the name of the tsv file
                    tsv_name=$(dx describe "$tsv_id" --json | jq -r '.name')

                    # compare the two files
                    echo "Comparing ${case_id}"

                    bash single_comparison.sh \
                        "$case_id" "$vcf_id" "$vcf_name" "$tsv_id" "$tsv_name" "$output_all_vars"
                else
                    echo "${case_id} not processed due to archived files"
                fi
            fi
        fi
    done

    i=$((i+1))
done

# run script to identify the cases each unique mismatch affects
echo "Identifying unique mismatch information"
python find_unique_vars.py "$output_all_vars" "$output_unique_vars"

# run script to apply panel filtering to all mismatches
echo "Identifying unique case information"
python apply_panel_filtering.py "$output_all_vars" "$output_panel_filtering"
