# vep_issue_testing

Confluence documentation: https://cuhbioinformatics.atlassian.net/wiki/spaces/RD/pages/3009019922/231009+Comparing+TSO500+CombinedVariantOutput.tsv+to+VEP-annotated+VCF

Objective:
For all Helios (TSO500) cases which annotated VCFs using VEP v107, compare the TSO500 output and the annotated VCF to identify any discrepancies in p. annotations; list all such discrepancies; and identify the number of cases which each unique discrepancy affects.

Usage:
- Executing unarchive_files.sh will identify all matched TSV and VCF files from relevant projects, and unarchive them if necessary.
- Executing process_all_cases.sh will identify all relevant projects in DNAnexus, then identify matched TSV and VCF files from that project.
- For each case, process_all_cases.sh will execute single_comparison.sh, which compares the relevant TSV and VCF files to identify variants with p. annotation discrepancies.
- This will generate a single output TSV file listing all instances of p. annotation discrepancies from every case.
- Executing find_unique_vars.py will condense this list, generating a list of unique annotation discrepancies along with the number and identity of the cases each affects.

Summary of results:
- 59 DNAnexus projects for Helios runs used VEP v107
- 779 unique cases for DNA analysis in these projects
- 694 cases had at least one variant where p. annotation differed between the TSO500 output and the annotated VCF
- 1843 individual instances of an annotation discrepancy across all cases
- 204 unique annotation discrepancies, affecting between 1 and 645 cases each
