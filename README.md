# TnClone
TnClone: High-throughput clonal analysis using Tn5-mediated library construction and _de novo_ assembly

TnClone now supports reference mapping based analysis

## NOTE
TnClone's manuscript is under preparation. Please do not copy/modify this code before publication.
## Before getting started
Please read throughly the manual of version user using. 
Almost all information including file formats for TnClone to run and output format description are described in manual.

## Getting started

### Prerequisites

JAVA 1.7 or higher<br>
Python 2.7<br>
matplotlib<br>
numpy<br>
HTSeq<br>

JAVA download and installation guides can be easily found at web.<br>

Please install PyQt4 following instructions here http://pyqt.sourceforge.net/Docs/PyQt4/installation.html<br>

To download python go to https://www.python.org/ and download and install python

TnClone was tested under Ubuntu 14.04 LTS and Windows 7

### Installation

#### Linux users

For Linux users there is no typical method to install this package.

Just clone this repository using command line below on your terminal screen.

```
git clone https://github.com/tahuh/tnclone.git
```
Or download this repository directly by clicking download button above and unzip the archive.

#### Windows users

Due to technical issue of Github page when uploading large data please visit the link below to download Windows version

Link : https://drive.google.com/drive/u/4/folders/10x1B7rVWPt0rPMyT2ZIMdAXSn-2VuhVu

Download all files and folders uploaded in the link above

### How to use

#### Test data
For assembly based method, use data under directory "denovo_assembly_data". 

For reference mapping module use data under "reference_mapping" directory

#### TnClone's different strategy for analysis
TnClone offers two different strategy for analysis currently one is de novo assembly based method and the other is reference mapping.

See the image attatched below

![ReferenceMapping](https://github.com/tahuh/tnclone/tree/master/image/fig2.png)

As shown in above image, for de novo assembly, user must check the assembly box in additional option panel.

If one wants to use reference mapping based method, user must uncheck the assembly box.

#### Linux users

After downloading TnClone please navigate to your directory where TnClone installed.

Here we asseume your path to TnClone is ~/Downloads/TnClone

Navigate to tnclone sub directory by typing

```
cd ~/Downloads/TnClone/tnclone/linux
```

if you have naviaged to the path above then type

```
python tnclone.py
```
Will give you a interface

User <b>MUST</b> execute above command where TnClone is installed!!!

For users who have experience using command line based analysis, we have tnclone-cmd.py under same directory where tnclone.py located

Below described options for TnClone's command line version - This will be deprecated soon. Use GUI version instead

```
Usage: tnclone-cmd.py [options]

TnClone CUI version

Options:
  -h, --help            show this help message and exit
  --ngs_path=NGS_PATH   Raw NGS result path. The original sequencing read file directory/path which is not sorted.
  --sort_file=SORT_FILE
                        Sorting information containing file.
  --sample_file=SAMPLE_FILE
                        Sample information file.
  --start=START         Start sequence of assembly. Length must be same as
                        k-mer size option
  --end=END             End sequence of assembly. Length must be same as k-mer
                        size option
  --ref_path=REF_PATH   Reference directory that contains reference files for
                        analysis
  --output_path=OUTPUT_PATH
                        Output path
  --kmer_size=KMER_SIZE
                        K-mer size. default = 63
  --min_kocc=MIN_KOCC   Minimum k-mer occurence.
                        This value is used to reduce ambiguity. default = 3
  --sort_on=SORT_ON     Check whether to turn on sorting. 0 for no 1 for yes.
                        default = 1
  --trim_on=TRIM_ON     Check whether to turn on trimming. 0 for no 1 for yes.
                        default = 1
  --assembly_on=ASSEMBLY_ON
                        Check whether to turn on assembly. 0 for no 1 for yes.
                        default = 1
  --analysis_on=ANALYSIS_ON
                        Check whether to turn on analysis. 0 for no 1 for yes.
                        default = 1
  --auto_analysis=AUTO_ANALYSIS
                        Generates analysis automatically. Will start from
                        analysis you specified
                        0 for no 1 for yes. default = 1
  --bed_file=BED_FILE   Bed file
  --down_sampling=DOWN_SAMPLING
                        Downsampling or not 1 for yes and 0 for not. 
                        default = 0
  --down_sample_rate=DOWN_SAMPLE_RATE
                        Down sampling ratio. default = 0.3
  --seed_mismatch_correction_method=SEED_MISMATCH_CORRECTION_METHOD
                        Type of methode to correct seed error. One of depth or
                        mismatch.
                        default=depth
  --num_mismatch=NUM_MISMATCH
                        Number of mismatch allowed.
                        if seed mismatch correction method is set depth , this
                        will not invoked.
                        default = 3
  --aln_match=ALN_MATCH
                        Match score. INT value. default = 1
  --aln_mismatch=ALN_MISMATCH
                        Mismatch score. INT value. Will work as negative
                        value. default = 3
  --aln_open_penalty=ALN_OPEN_PENALTY
                        Gap open penalty. INT value. default = 5
  --aln_gap_extension=ALN_GAP_EXTENSION
                        Gap extension penalty. INT value. Default = 3
  --draw_coverage=DRAW_COVERAGE
                        Draw coverage plot before analyis. 0 for no 1 for yes.
                        default = 0
```



After downloading TnClone one <b>MUST</b> go to your directory where TnClone is installed

Here we asseme your directory is  C:\\User\\Downloads

Search for folder C:\\User\\Downloads\\windows\\tnclone\\dist

There is an executable file called 'tnclone'.

By double clicking this icon, tnclone will execute with some command line prompt automatically open.

## Output format
Basically TnClone generates assembled contigs in "assem" directory and file name contain "final" is final assembly sequence.
If TnClone has reference sequence as input TnClone will generate summary file in "report" directory and variant information in "vcf" directory.
More detailed aligned information are located under "sam" directory for user to use if reference sequence is provided.

Summary output format is described below
```
[Contig Number Summary]
Q03JI6	1(C=1:I/D=0)	1(C=1:I/D=0)	0(C=0:I/D=0)	0(C=0:I/D=0)
[DNA ERROR/FREE SUMMARY]
Q03JI6	(E)	(E)	(NOINFO)	(NOINFO)
[PROTEIN ERROR/FREE SUMMARY]
Q03JI6	(E/E/T)	(E/E/T)	(NOINFO)	(NOINFO)
```
The output summary consist of three fields (contig number summary, DNA error free summary, protein error free summary).
For contig number summary the first column is always the gene name specified for TnClone's input and the number of samples analyzed are follows.
For each column following gene column, the first number is total number of contigs of assembly. The numbers in parenthesis is detailed description of contigs. Letter with "C" is "confident sequence", letter with "I/D" is "indel sequence" and letter with "NC" is "not confident sequence". 
For dna/protein error free summary field, character "E" tells given contig/protein is error compared to reference and "T" tells given contig/protein is error free.

Variant format is described below
```
#VCF(simple version.)
#INFO:The 'short' mark in INFO field means the query is shorter than reference. INDEL or large soft clipping
#INFO:The 'homopoly' mark in INFO field means there is homopolymer region in refrence. These variant might be caused by sequencing errors.
#REF	QUERY	POS	REF	ALT	INFO
A0Q5Y3	contig:375	1600	A	C	SHORT
A0Q5Y3	contig:375	2368	T	C	SHORT
A0Q5Y3	contig:375	2731	T	C	SHORT
A0Q5Y3	contig:375	3528	C	T	SHORT
A0Q5Y3	contig:375	4269	T	C	SHORT
A0Q5Y3	contig:375	4848	A	C	SHORT
A0Q5Y3	contig:375	4879	C	T	SHORT
A0Q5Y3	contig:141	2368	T	C	SHORT
A0Q5Y3	contig:141	2731	T	C	SHORT
A0Q5Y3	contig:141	3528	C	T	SHORT
A0Q5Y3	contig:141	4269	T	C	SHORT
A0Q5Y3	contig:141	4848	A	C	SHORT
A0Q5Y3	contig:141	4879	C	T	SHORT
```
TnClone mimics VCF format as variant output but modified way. Each line consist of reference sequence name, contig name, variant position, reference base, contig base and extra information.
If indel present, reference will be dash(-) for insertion and contig base will be dash(-) for deletion.

See more detailed information in manual text.
## Test data
Use test data in data directory.
FASTQ files are already trimmed so user only need to check assembly and sort option and turn off sort and trim option

## Some troubleshoot solution
### TnClone FAILs to assembly!!!
If TnClone fails to assembly, please check if your sequencing depth is enough. <br>
In TnClone paper, assembly was failed when sequening depth is less than 400(< 400X). <br><br>
To solve this problem here are two options one can try <br>
1. Use Diagnosis plot to see if (1) your sequencing depth was enough (2) the region of interest is properly sequenced.<br>
2. Change size of k-mer smaller than current one(such as 55 or 47. Use odd number)<br>
3. If 2 fails then change min kmer occurence in Additional options tab to zero<br>
### TnClone shows NO SEED in output file!!!
If one see NO SEED in output(usually 8 bytes file size of .contig file), then please select other region of seeds <br>
Please avoid repetitive region as seed/terminal sequence<br>
### TnClone result has 0 bytes outfile
This isssue says your assembly has failed. Please see <b>TnClone FAILs to assembly</b> section</br>
### TnClone cannot recognize the path(directory) for required fields
There are two main reasons that TnClone cannot recognize the path(directory)<br>
1. There is *Really NO such path*
2. Please check if you have typed the path directly to the screen.<br>
   TnClone ONLY accepts required directories selected by *BROWSE* button<br>
3. For windows users please check if your path separator ("\\" sign) is set twice(i.e. \\\\).<br>
### TnClone says "QThread: Program exited while thread is running." when exit the program
This issue seems okay when the developers analysed data but will be fixed soon.
### Diagnosis does not work!
TnClone utilized bwa and samtools for diagnosis which is installed in the exact same path where tnclone.py lies.<br>
So please *Execute TnClone where tnclone.py is installed*<br>
### Trmming does not work!
As same as above please *execute TnClone where tnclone.py is installed*
### Even I have executed TnClone where the tnclone.py, JNI Error is occurred or diagnosis/trimming do not occur!
This is usually happens on Linux system since user do not have permission to execute auxillary files required TnClone to perform.<br>
Please type the command where the TnClone is installed<br>
As mentioned above, forget about typing the '$' sign
```
/Directory/Where/You/Have/Installed/TnClone$ chmod +x wakeup.sh
/Directory/Where/You/Have/Installed/TnClone$ ./wakeup.sh
```
Eventhough you have tested above procedure and gets error, please check https://mafft.cbrc.jp/alignment/software/installation_without_root.html
for complie mafft<br>
Also, if JNI error occurs then navigate msa2vcf at 
http://lindenb.github.io/jvarkit/
Then download jvarkit and follow installation instructions there.<br>
Copy msa2vcf.jar by command line below<br>
Suppose you have downloaded jvarkit at /home/username/Downaloads and TnClone is under /home/username/Downloads/tnclone
```
$ cp /home/username/Downloads/dist/msa2vcf.jar /home/username/Downloads/tnclone/linux
```

Then perform TnClone<br>
Hope this helps!!!!
## Authors
* Sunghoon Heo - *Initially developed/maintaining this software* -
* Byungjin Hwang  - *Contributor* -

## License
TnClone - GPL v3
GATK - BSD 3-clause "New" or "Revised" license
Picard - MIT

## Acknowledgements
We thanks to Bang's lab for kind donation of their samples for testing TnClone.

We also thanks to Junho Jung's lab for kind donation of their scFv samples.

## Contact
For any other issue while using TnClone, first open issue of TnClone. 
If issue is not helpful please contact via team.tnclone2@gmail.com
