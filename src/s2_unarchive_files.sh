#!/bin/bash

762 cases, each with two files (tsv and vcf)


## IDENTIFY ARCHIVED FILES

while read -r case_id vcf_name vcf_id tsv_name tsv_id; do
    for file_id in "$vcf_id" "$tsv_id"; do

        state=$(dx describe "$file_id" --json | jq -r '.archivalState')

        if [[ "$state" == "archived" ]]; then
            echo "$file_id" >> files/archived_files

        fi
    done
done < files/s1_f6_final_output


## UNARCHIVE ARCHIVED FILES

while read -r file_id; do

    dx unarchive "$file_id" -y

done < files/archived_files


## CHECK FILE UNARCHIVING

while read -r file_id; do

    state=$(dx describe "$file_id" --json | jq -r '.archivalState')
    echo "$state"

done < files/archived_files
