# CNV-PCC

CNV-PCC: Detection of copy number variations from Next-generation Sequencing Data

1. Installation:

Basic requirements:

Software: Python, R, SAMtools, BWA

Operating System: Linux

Python version: 3.8.5 and the higer version

R version: 4.0.4 and the higer version

Required python packages: 

- numpy
-  pandas
-  pysam
- subprocess
-  pyod
-  numba
- imblearn

Required R packages:

- DNAcopy

2. Running software:

2.1 Preprocessing of input files:

Usually, the following documents are required:

A genome reference sequence fasta file.
A bam file from a  sample.

The bam file  must be indexed. You can do the following:
$samtools index example_sorted.bam

2.2 Operating command:

python CNV-PCC.py [bamfile] [reference] [binsize] [Ls] [outfile] 

bamfile: The path to the bam file representing the sample used by the user.

reference: The path to the fasta file of the genome reference sequence used by the user.

binsize: The size of bin. The recommended value is 1000.

Ls: The length of subsegment. The recommended value is 10000.

outfile: The path to the output file, which contains the results of the detected CNVs.



