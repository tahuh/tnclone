#!/usr/bin/python

"""
Trimmer.py

Trimming Illumina reads using Trimmomatic

Author Sunghoon H.

"""
import subprocess
import commands
import os
from SequenceParser import FQParser
class Trimmer(object) :
	def __init__(self, idir, odir, conf, trimmo_path,is_cmd = False) :
		self.idir = idir
		self.odir = odir
		self.conf = conf
		self.tpath = trimmo_path
		#self.editor = editor
		if is_cmd :
			self.path = os.path.dirname(os.path.abspath(__file__)) + "/Trimmomatic-0.36/adapters/"
		else:
			self.path = ''
			
	def _load_fq(self, parser):
		d = {}
		for id , seq, qual in parser.parse():
			d[id.split()[0]] = (id, seq, qual)
		return d
		
	def _build_fq_block(self, id , data):
		id , seq , qual = data[id]
		return id + "\n" + seq + "\n+\n" + qual
	def trim(self) :
		input_dir = self.idir
		output_dir = self.odir
		trimopath = self.tpath
		with open(self.conf) as config_file:
			for line in config_file :
				data = line.rstrip().split("\t")
				gene = data[0]
				if not os.path.isdir(output_dir + "/" + gene) :
					os.mkdir(output_dir + "/" + gene)
				n_sample = int(data[1])
				for i in xrange(1,n_sample+1) :
					trimmoscript = "java -jar " + trimopath + " PE -phred33 " + \
						       input_dir + "/" + gene + "/" + gene + "_" + str(i) + "_R1.sorted.fastq " + \
						       input_dir + "/" + gene + "/" + gene + "_" + str(i) + "_R2.sorted.fastq " + \
						       output_dir + "/" + gene + "/" + gene + "_" + str(i) + "_R1.trimmed.fastq.tmp " + \
						       output_dir + "/" + gene + "/" + gene + "_" + str(i) + "_R1_unpaired.fq.tmp.gz " + \
						       output_dir + "/" + gene + "/" + gene + "_" + str(i) + "_R2.trimmed.fastq.tmp " + \
						       output_dir + "/" + gene + "/" + gene + "_" + str(i) + "_R2_unpaired.fq.tmp.gz " + \
						       "ILLUMINACLIP:%sTruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"%(self.path)
					#print trimmoscript
					#print
					#self.editor.append(trimmoscript + "\n")
					subprocess.call(trimmoscript,shell=True)
					#commands.getoutput(trimmoscript)
					
					
					### Trim out N containing reads
				
					#self.editor.append("[TnClone::Trimmer] TnClone is trimming additional N containing reads.")
					#self.editor.append("[TnClone::Trimmer] Memory might be inflating due to internal memory consumption.")
					fqparser1 = FQParser(output_dir + "/" + gene + "/" + gene + "_" + str(i) + "_R1.trimmed.fastq.tmp")
					fqparser2 = FQParser(output_dir + "/" + gene + "/" + gene + "_" + str(i) + "_R2.trimmed.fastq.tmp")
					fqparser1.open()
					fqparser2.open()
					fq1_data = self._load_fq(fqparser1)
					fq2_data = self._load_fq(fqparser2)
					
					fqparser1.close()
					fqparser2.close()
					id_set = fq1_data.keys()
					
					discards = []
					total_reads = len(id_set)
					for id in id_set:
						seq1 = fq1_data[id][1]
						seq2 = fq2_data[id][1]
						
						if ('N' in seq1) or ('N' in seq2) :
							discards.append(id)
							
					for id in discards :
						fq1_data.pop(id)
						fq2_data.pop(id)
						
					
					sorted_keys = fq1_data.keys()
					
					sorted_keys.sort()
					remain_set = len(sorted_keys)
					otrimmed1  =  open( output_dir + "/" + gene + "/" + gene + "_" + str(i) + "_R1.trimmed.fastq", "w")
					otrimmed2  =  open( output_dir + "/" + gene + "/" + gene + "_" + str(i) + "_R2.trimmed.fastq", "w")
					for key in sorted_keys:
						fq1_block = self._build_fq_block(key ,fq1_data)
						fq2_block = self._build_fq_block(key ,fq2_data)
						otrimmed1.write(fq1_block + "\n")
						otrimmed2.write(fq2_block + "\n")
						
					otrimmed1.close()
					otrimmed2.close()
					#try:
					#	self.editor.append("[TnCloneTrimmer] From %d read pairs %d read pairs are removed (%8.3f %%)"%(total_reads, remain_set , float(remain_set) / total_reads * 100))
					#except ZeroDivisionError:
					#	self.editor.append("[TnCloneTrimmer] No reads detected.")
