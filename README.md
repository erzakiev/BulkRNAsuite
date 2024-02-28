# BulkRNAsuite

An imagination of a start-to-finish suite for Bulk RNAseq analysis that 
- starts with user-submitted FASTQ files (uploaded by the user)
- checks the consistency of filenames
   - finds and merges multiplexed files, if applicable (_001, _002, _003, _004)
   - finds and treats together R1 and R2 files, if applicable
- runs Fastqc
- runs bbduk.sh to trim adapters
- runs Fastqc again on the trimmed files
- aligns the trimmed files with Salmon using GRCm39 or GRCh38 genome full indices with decoys (SAF)

## Reporting
- QC results before trimming
- QC results after trimming
- TPM matrix
- DEGs for each pairwise comparison in the list of all available groups
- PCA
- Interactive (plotly-based) Volcano plots for each pairwise comparison in the list of all available groups
- Enriched terms based on the DEGs for each pairwise comparison in the list of all available groups
- Dorothea TF regulon enrichment for each pairwise DEGs
- GSVA analysis
