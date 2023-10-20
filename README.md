# vep_issue_testing

Confluence documentation: https://cuhbioinformatics.atlassian.net/wiki/spaces/RD/pages/3009019922/231009+Comparing+TSO500+CombinedVariantOutput.tsv+to+VEP-annotated+VCF

These scripts were written in order to:
- Identify all Helios (TSO500) cases which involved DNA analysis using VEP v107
- List the CombinedVariantOutput.tsv file (the raw output of the TSO500 app) and the VEP-annotated VCF for each of these cases
- For each case, compare the TSO500 output and the annotated VCF to identify any discrepancies in p. annotations
- List all such discrepancies
- Identify the number of cases which each unique discrepancy affects

The purpose of each script is as follows:
- s1_get_cases_and_files.sh (step 1): Generate text files listing the cases involved, and their output files to compare
- s2_unarchive_files.sh (step 2): For each of the output files identified in step 1, identify whether or not the file needs unarchiving in DNAnexus
- s3_compare_files_single.sh (step 3): For a single case, compare the tsv and vcf output files and return a list of variants with annotation discrepancies. Used for testing the comparison process before applying it to all cases.
- s4_compare_files_all.sh (step 4): Apply the tsv/vcf comparison process to all cases listed as part of step 1, thereby identifying all instances of variant annotation discrepancies across all cases
- s5_find_unique_vars.py (step 5): Condense the list of annotation discrepancies generated in step 4, generating a list of unique annotation discrepancies along with the number and identity of the cases each affects

Key outputs:
- files/s1_f1_tso_projects: Lists DNAnexus projects for Helios runs which used VEP v107
- files/s1_f5_dna_case_list: Lists all unique samples from those cases for which DNA analysis was performed
- files/s1_f6_final_output: Lists the TSV and VCF file names and IDs for each of the DNA analysis samples
- files/cvo_vs_vep_output_all: Lists all instances of annotation discrepancies between TSV and VCF output for all cases
- files/unique_mismatch_variants: Lists all unique annotation discrepancies, and the number and identity of the cases they affect

Summary of results:
- 59 DNAnexus projects for Helios runs which used VEP v107
- 779 unique cases for DNA analysis in these projects
- 694 cases had at least one variant where p. annotation differed between the TSO500 output and the annotated VCF
- 1843 individual instances of an annotation discrepancy across all cases
- 204 unique annotation discrepancies, affecting between 1 and 645 cases each
