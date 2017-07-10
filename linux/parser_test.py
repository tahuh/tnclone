#!/usr/bin/python 

import os ### for pid
import time ### for time check
import gzip ### for benchmark with biopython

from Bio import SeqIO

from SequenceParser import FAParser
from SequenceParser import FQParser

from memory import MemoryUsageManager


### FASTA
faname = "/Data/reference/reference/hg19.fa"

### FQ in gzip compressed
gz_fqname = "/Data/HSH/projects/cas9/ALL_PLATFORM_RERUN_160810/160708_HiSeq/raw/TN1607L0001--AGCGATAG-GAACCATT_S110_L008_R1_001.fastq.gz"

### FQ uncompressed
fqname = "/Data/HSH/projects/cas9/ALL_PLATFORM_RERUN_160810/160708_HiSeq/trimmed/A0Q5Y3/A0Q5Y3_6_R1.trimmed.fastq"



pid = os.getpid()
manager = MemoryUsageManager(pid, unit="Gb")

############################ FAPARSER TEST ##################################
fa_dict = {}
print "Step1 Using method 1 of synopsis for fasta"
faparser = FAParser(faname)

faparser.open()

for id , desc , seq in faparser.parse() :
	#print id
	fa_dict[id] = seq

faparser.close()

del faparser
manager.profile()

print manager.mem_usage

save1 = manager.save_point
fa_dict.clear()

print "Step2 Using method 2 of synopsis for fasta"

faparser = FAParser()

faparser.open(faname)

t1 = time.time()
for id , desc , seq in faparser.parse():
	#print id
	fa_dict[id] = seq
t2 = time.time()

diff_custom = t2 - t1

faparser.close()
del faparser

manager.profile()

print manager.mem_usage
save2 = manager.save_point

print "Mem diff = %f Gb ( expected 0 ) " %(save2-save1)

fa_dict.clear()

with open(faname) as fasta_file:
	t1 = time.time()	
	for record in SeqIO.parse(fasta_file, "fasta") :
		id = record.id
		seq = str(record.seq)
		fa_dict[id] = seq
	t2 = time.time()
	seqio_time = t2-t1
manager.profile()
save3 = manager.save_point

print "Custom Parser Elapsed                    SeqIO Elapsed"
print "%f sec                                   %f sec"%(diff_custom , seqio_time)

print "Custom parser Memusage                   SeqIO Mem usage"
print "%f Gb                                    %f Gb"%(save2,save3)

fa_dict.clear()


######################### FQPARSER TEST #################################

fq_dict = {}

print "Step 3 Method 1 in the synopsis for fastq ( gzip file )"

fqparser = FQParser()

fqparser.open(gz_fqname)

fqparser.open()

t1 = time.time()
for id , seq, qual in fqparser.parse():
	fq_dict[id] = {"seq" : seq , "qual": qual}
t2 = time.time()
fqparser.close()

manager.profile()

save = manager.save_point
#print manager.mem_usage


custom_gz_fq_time = t2 - t1

fq_dict.clear()

print "Compare same job with SeqIO module"

with gzip.open(gz_fqname) as gz_fqfile:
	t1 = time.time()
	for record in SeqIO.parse(gz_fqfile , "fastq") :
		id = record.id
		seq= str(record.seq)
		qual = "".join([chr(x+33) for x in record.letter_annotations["phred_quality"]])
		fq_dict[id] = {"seq" : seq, "qual" : qual}
	t2 = time.time()

seqio_gz_fq_time = t2-t1
manager.profile()
save2 = manager.save_point

fq_dict.clear()

print "Mem usage by Custom Parser           Mem usage by SeqIO"
print "%f Gb                                %f Gb"%(save,save2)
print "Time elapsed by custom parser        Time elapsed by SeqIO"
print "%f sec                               %f sec"%(custom_gz_fq_time,seqio_gz_fq_time)

print "Step 4 Method 2 in the synopsis for fastq ( gzip file )"

fqparser = FQParser(gz_fqname)

fqparser.open()

fqparser.open()

t1 = time.time()
for id , seq, qual in fqparser.parse():
	fq_dict[id] = {"seq" : seq , "qual": qual}
t2 = time.time()
fqparser.close()

manager.profile()

save = manager.save_point
#print manager.mem_usage


custom_gz_fq_time = t2 - t1

fq_dict.clear()

print "Compare same job with SeqIO module"

with gzip.open(gz_fqname) as gz_fqfile:
	t1 = time.time()
	for record in SeqIO.parse(gz_fqfile , "fastq") :
		id = record.id
		seq= str(record.seq)
		qual = "".join([chr(x+33) for x in record.letter_annotations["phred_quality"]])
		fq_dict[id] = {"seq" : seq, "qual" : qual}
	t2 = time.time()

seqio_gz_fq_time = t2-t1
manager.profile()
save2 = manager.save_point
fq_dict.clear()
print "Mem usage by Custom Parser           Mem usage by SeqIO"
print "%f Gb                                %f Gb"%(save,save2)
print "Time elapsed by custom parser        Time elapsed by SeqIO"
print "%f sec                               %f sec"%(custom_gz_fq_time,seqio_gz_fq_time)

print "Step 5 Method 1 in the synopsis for fastq ( uncompressed file )"

fqparser = FQParser()

fqparser.open(fqname)

fqparser.open()

t1 = time.time()
for id , seq, qual in fqparser.parse():
	fq_dict[id] = {"seq" : seq , "qual": qual}
t2 = time.time()
fqparser.close()

manager.profile()

save = manager.save_point
#print manager.mem_usage


custom_fq_time = t2 - t1

fq_dict.clear()

print "Compare same job with SeqIO module"

with open(fqname) as fqfile:
	t1 = time.time()
	for record in SeqIO.parse(fqfile , "fastq") :
		id = record.id
		seq= str(record.seq)
		qual = "".join([chr(x+33) for x in record.letter_annotations["phred_quality"]])
		fq_dict[id] = {"seq" : seq, "qual" : qual}
	t2 = time.time()

seqio_fq_time = t2-t1
manager.profile()
save2 = manager.save_point

fq_dict.clear()

print "Mem usage by Custom Parser           Mem usage by SeqIO"
print "%f Gb                                %f Gb"%(save,save2)
print "Time elapsed by custom parser        Time elapsed by SeqIO"
print "%f sec                               %f sec"%(custom_fq_time,seqio_fq_time)

print "Step 6 Method 2 in the synopsis for fastq ( uncompressed file )"

fqparser = FQParser(fqname)

fqparser.open()

fqparser.open()

t1 = time.time()
for id , seq, qual in fqparser.parse():
	fq_dict[id] = {"seq" : seq , "qual": qual}
t2 = time.time()
fqparser.close()

manager.profile()

save = manager.save_point
#print manager.mem_usage


custom_fq_time = t2 - t1

fq_dict.clear()

print "Compare same job with SeqIO module"

with open(fqname) as fqfile:
	t1 = time.time()
	for record in SeqIO.parse(fqfile , "fastq") :
		id = record.id
		seq= str(record.seq)
		qual = "".join([chr(x+33) for x in record.letter_annotations["phred_quality"]])
		fq_dict[id] = {"seq" : seq, "qual" : qual}
	t2 = time.time()

seqio_fq_time = t2-t1
manager.profile()
save2 = manager.save_point
fq_dict.clear()

print "Mem usage by Custom Parser           Mem usage by SeqIO"
print "%f Gb                                %f Gb"%(save,save2)
print "Time elapsed by custom parser        Time elapsed by SeqIO"
print "%f sec                               %f sec"%(custom_fq_time,seqio_fq_time)
              
