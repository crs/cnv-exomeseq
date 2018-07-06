
### CNV from exomeseq Snakemake pipeline

usage:
	snakemake --snakefile Snakefile.py

Req.:
- samtools
- CNVnator
- Python

* v0.2
	- Changed pipeline to work on fastq.gz
	- included binsize parameter to output files

* v0.1 
	- Single sample pipeline, folder tracking not yet implemented but intended.
	- Works on paired-end data (R1, R2) and two lanes for each sample
