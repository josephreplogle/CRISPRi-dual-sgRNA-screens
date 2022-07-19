# CRISPRi-dual-sgRNA-screens

This repository contains scripts for alignment of sequencing data from dual-sgRNA CRISPR screening data. These scripts were adapted from: https://github.com/mhorlbeck/ScreenProcessing

For alignment of data without UMI, use python dualguide_fastqgz_to_counts.py. Example command:

```bash
python dualguide_fastqgz_to_counts.py 20200513_library_1_2_unbalanced.csv  alignments fastq/UDP0007* fastq/UDP0006*
```

Here, fastq/UDP0007* and fastq/UDP0006* are fastq files (the script will automatically detect read 1 and read 2; you can add as many fastq files as you want); "alignments" is the output directory; 20200513_library_1_2_unbalanced.csv is the library file. The output counts file can then be input into MAGECK or analyzed directly.


For libraries with UMI, demultiplexing on only the i5 index using the i7 index (IBC) as a read is performed as detailed: https://gist.github.com/sumeetg23/a064a36801d2763e94da2e191699fb9f. For alignment of data with UMI, use python dualguide_UMI_fastqgz_to_counts.py. Example command: 

```bash
python dualguide_UMI_fastqgz_to_counts.py 20200513_library_1_2_unbalanced.csv \
  UMI_sequences_filtered.csv \
  alignments \
  JR72_S1_L001* JR72_S1_L002* JR72_S1_L003* JR72_S1_L004*
```
This will output a file of counts for each sgRNA-UMI combination. From this, we can then calculate the enrichment of each guide between screen arms in python or MaGECK.



