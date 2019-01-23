# TnClone
TnClone - High-throughput clonal analysis using Tn5-mediated library construction and de novo assembly

# How to test data

Create an output directory such as output

Then create subdirectory called "trimmed"

Below trimmed create two directories "A8REH9" and "A5EIM8" and locate fastq files into those directories

Namely, locate "A8REH9_1_R1.trimmed.fastq" and "A8REH9_1_R2.trimmed.fastq" into "A8REH9" directory, and perform other files likely in this manner.

Create a directory to store reference sequence called "ref"

Put all reference sequence related files (.fa , .fa.dict, .fa.xxx) in this directory

Use "region.bed" to create assembly target for TnClone to assemble (used in additional option tab. see manual).

Then execute with de novo assembly engine workflow as decribed in manual.
