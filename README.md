# RNAgenotyping_hybrids
Scripts for genotpying analysis using RNAseq data of F1 hybrids of the grass Alloteropsis semialata

# Pipeline for SNP genotyping using RNAseq data, designed to determine the pollen parents of hybrids between photosynthetic types in Alloteropsis semialata 

## Part 1. Generate genotypes file (see Methods section in the paper for details)
# 1. Map reads to genome of Alloteropsis semialata (Dunning et al. 2019; CDS file) 
# 2. Generate VCF file after filtering the bam alignment
# 3. Filter VCF
# 4. Retain only multiallelic sites in which all hybrids are heterozygous and all potential parents are homozygous
# 5. Format the table before R analysis -> Script: "1_format_SNP_table.sh"; output -> "variants_only_multiallelic_sites_heter_and_homoz_reformatted.txt"

## Part 2. Genotyping analysis using R
# Input 1: "variants_only_multiallelic_sites_heter_and_homoz_reformatted.txt" # generated in Part 1
# Input 2: "metadata.csv" # to get sample info

# Script: "2_genotyping_hybrids_with_rna.R"
