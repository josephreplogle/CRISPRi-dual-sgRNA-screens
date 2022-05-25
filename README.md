# CRISPRi-dual-sgRNA-screens

This repository contains scripts for alignment of sequencing data from dual-sgRNA CRISPR screening data. These scripts were adapted from: https://github.com/mhorlbeck/ScreenProcessing

For alignment of data without UMI, use python dualguide_fastqgz_to_counts.py. Example command:

python dualguide_fastqgz_to_counts.py 20200513_library_1_2_unbalanced.csv  alignments fastq/UDP0007* fastq/UDP0006*

Here, fastq/UDP0007* and fastq/UDP0006* are fastq files (the script will automatically detect read 1 and read 2; you can add as many fastq files as you want); "alignments" is the output directory; 20200513_library_1_2_unbalanced.csv is the library file. The output counts file can then be input into MAGECK or analyzed directly.


For alignment of data without UMI, use python dualguide_UMI_fastqgz_to_counts.py. Example command: 

python dualguide_UMI_fastqgz_to_counts.py 20200513_library_1_2_unbalanced.csv \
  UMI_sequences_filtered.csv \
  alignments \
  JR72_S1_L001* JR72_S1_L002* JR72_S1_L003* JR72_S1_L004*

This will output a file of counts for each sgRNA-UMI combination. From this, we can then calculate the enrichment of each guide between screen arms in python or MaGECK.



