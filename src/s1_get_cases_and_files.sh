#!/bin/bash


file_1="files/s1_f1_tso_projects"
file_2="files/s1_f2_dx_info"
file_3="files/s1_f3_all_cases_and_files"
file_4="files/s1_f4_dna_cases_only"
file_5="files/s1_f5_dna_case_list"
file_6="files/s1_f6_final_output"

## LIST ALL AFFECTED PROJECTS (59 projects)

printf "\nListing affected projects\n\n"

dx find projects --name "002_23052*_TSO500" > "$file_1"  # May
dx find projects --name "002_2306*_TSO500" >> "$file_1"  # June
dx find projects --name "002_2307*_TSO500" >> "$file_1"  # July
dx find projects --name "002_2308*_TSO500" >> "$file_1"  # Aug
dx find projects --name "002_2309*_TSO500" >> "$file_1"  # Sept


## LIST FILE DATA (820 tso500 output tsvs, 763 annotated vcfs total)

printf "\nListing files in affected projects\n\n"

while read -r project; do

    for suffix in \
    *withLowSupportHotspots_annotated.vcf.split.vcf.gz \
    *CombinedVariantOutput.tsv; do

        pid=$(echo "$project" | cut -f 1 -d " ")

        dx find data \
        --name "${suffix}" \
        --project "$pid" \
        --json | jq -r '.[] | .project + ":" + .id + "\t" + .describe.name'\
        >> "$file_2"

    done
done < "$file_1"


## GET CASE IDS (820 unique cases, 758 for DNA analysis)

printf "\nGetting case ids and removing duplicates\n\n"

while read -r file_id file_name; do

    case_prefix=$(echo "$file_name" | cut -f 1 -d $'_')

    # exclude cases with Epic RNA analysis ids
    if [[ $file_name != *"-8472"* ]] && [[ $file_name != *"-9689"* ]]; then

        case="$case_prefix"

        # some older cases have awkward prefixes
        if [[  $case_prefix == *"egg6" ]]; then
            case=$(echo "$file_name" | cut -f 1 -d $'-')
        fi

        printf "%s\t%s\t%s\n" "$case" "$file_name" "$file_id" \
        >> "$file_3"

    fi
done < "$file_2"

# remove duplicate records, sort by case id, get list of unique dna cases
awk '!seen[$2]++' "$file_3" | sort -u -k1 > "$file_4"
awk -F"\t" '{print $1}' "$file_4" | uniq > "$file_5"


## LIST THE 2 FILES FOR EACH UNIQUE DNA CASE

printf "\nCollating file data for each case\n\n"

while read -r case_id; do

    vep_vcf_name=$(awk -F"\t" -v case_num="$case_id" \
    '$1==case_num && $2~"vcf.gz" {print $2}' "$file_4")

    vep_vcf_id=$(awk -F"\t" -v case_num="$case_id" \
    '$1==case_num && $2~"vcf.gz" {print $3}' "$file_4")

    tso_output_name=$(awk -F"\t" -v case_num="$case_id" \
    '$1==case_num && $2~".tsv" {print $2}' "$file_4")

    tso_output_id=$(awk -F"\t" -v case_num="$case_id" \
    '$1==case_num && $2~".tsv" {print $3}' "$file_4")

    if [[ "$vep_vcf_name" != "" ]]; then
        printf "%s\t%s\t%s\t%s\t%s\n" "${case_id}" "${vep_vcf_name}" \
        "${vep_vcf_id}" "${tso_output_name}" "${tso_output_id}" \
        >> "$file_6"
    fi

done < "$file_5"
