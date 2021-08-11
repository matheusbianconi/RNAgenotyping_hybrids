<h1> Genotyping F1 hybrids using RNAseq data </h1>
Scripts for genotpying analysis using RNAseq data of F1 hybrids of the grass <i>Alloteropsis semialata</i> from 

<b>Bianconi ME, Sotelo G, Curran EV, Milenkovic V, Samaritani E, Dunning LT, Osborne CP, Christin P-A.</b> 
Upregulation of C<sub>4</sub> characteristics does not consistently improve photosynthetic performance in intraspecific hybrids of a grass. 
<i>bioRxiv 2021.08.10.455822</i>; doi: https://doi.org/10.1101/2021.08.10.455822 

<h4> Part 1. Generate genotypes file (see Methods section in the manuscript for details) </h4>

<b>1)</b> Map filtered FASTQ reads to the genome of <i>Alloteropsis semialata</i> (CDS file) (Dunning et al. 2019)

<b>2)</b> Generate VCF file and filter it

<b>3)</b> Retain only multiallelic sites in which all F1 hybrids are heterozygous and all potential parents are homozygous

<b>4)</b> Reformat the VCF file before R analysis:
Use the script <code>1_format_SNP_table.sh</code> to generate the file <code>variants_only_multiallelic_sites_heter_and_homoz_reformatted.txt</code>

<h4> Part 2. Genotyping analysis using R </h4>

Use the R script: <code>2_genotyping_hybrids_with_rna.R</code>

Input files:

<b>1)</b> <code>variants_only_multiallelic_sites_heter_and_homoz_reformatted.txt</code> (generated above)

<b>2)</b> <code>metadata.csv</code> (sample info)


