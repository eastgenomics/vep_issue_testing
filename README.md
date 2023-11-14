# vep_issue_testing

Confluence documentation: https://cuhbioinformatics.atlassian.net/wiki/spaces/RD/pages/3009019922/231009+Comparing+TSO500+CombinedVariantOutput.tsv+to+VEP-annotated+VCF

## Objective

For all Helios (TSO500) cases which annotated VCFs using VEP v107:
- Compare the TSO500 output and the annotated VCF to identify any discrepancies in p. annotations
- List all instances of such discrepancies
- Identify the number of cases that each unique discrepancy affects
- By comparing affected genes to case panels, identify whether any mismatches could have affected diagnosis

## Usage

Executing process_all_cases.sh will:

1. Identify all relevant projects in DNAnexus, then identify matched TSV and VCF files from each project.
2. Execute single_comparison.sh once per case. This script compares the TSV and VCF files, and identifies variants with p. annotation discrepancies.
3. Generate an output file (output_all_mismatches) listing all instances of p. annotation discrepancies from every case.
4. Execute find_unique_vars.py. This script condenses output_all_mismatches, creating a list of unique mismatches along with the cases they each affect (output_unique_mismatches).
5. Execute apply_panel_filtering.py. This script uses two mapping files ('case_phenotypes.tsv' and 'gene_panels_map.txt') to identify whether each mismatch instance affects a gene which is actually in the panel used for that case (output_panel_filtering). The phenotype (clinical indication) for each case was identified from Epic.

There are 70 DNAnexus projects of the form 002_*_TSO500, each with up to 16 DNA samples, so there are a large number of potentially affected cases to process. This means that this code has a long run time (4-5 hours).
