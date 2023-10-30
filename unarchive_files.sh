#!/bin/bash

# identify all 002_*_TSO500 projects between 22 May and 20 Oct 2023
projects=$(dx find projects \
--name "002_*_TSO500" \
--created-after 2023-05-20 \
--created-before 2023-10-20 \
--brief)

for project in $projects; do

    # find IDs of VEP-annotated VCF files
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

        # identify the TSO500 output file
        # same tsv file is in project twice in different folders, only keep one
        tsv_id=$(dx find data \
            --name "${case_id}_CombinedVariantOutput.tsv" \
            --project "$project" \
            --brief | cut -d $'\n' -f 1)

        # if there is a matching tsv file,
        if [[ "${tsv_id// /}" != "" ]]; then

            tsv_name=$(dx describe "$vcf_id" --json | jq -r '.name')

            # check whether the vcf or tsv are archived
            vcf_state=$(dx describe "$vcf_id" --json | jq -r '.archivalState')
            tsv_state=$(dx describe "$tsv_id" --json | jq -r '.archivalState')

            # unarchive files if needed
            if [[ "$tsv_state" == 'archived' ]]; then
                echo "${tsv_name} is archived, unarchiving"
                dx unarchive "$tsv_id"
                sleep 2
            fi
            if [[ "$vcf_state" == 'archived' ]]; then
                echo "${vcf_name} is archived, unarchiving"
                dx unarchive "$vcf_id"
                sleep 2
            fi
        fi
    done
done
