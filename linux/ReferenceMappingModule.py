#!/usr/bin/python

"""
ReferenceMappingModule.py

Author : Thomas Sunghoon Heo

Using raw reads to map on reference and analyse data
"""
import os
import subprocess
class RawReadMappingEngine:
	def __init__(self, sample_name, fq1name, fq2name, refname, samname):
		self.sample_name = sample_name
		self.fq1name = fq1name
		self.fq2name = fq2name
		self.refname = refname
		self.samname = samname
		self.bamname = self.samname + ".bam"
		self.sortbam = self.samname + ".bam.sort"
		self.sortbam_rgadd = self.samname + ".bam.sort.bam.rgadd.bam"
	def runMapping(self):
		indexing =  "./bwa index %s"%(self.refname)
		samtools_faidx = "./samtools faidx %s"%(self.refname)
		command = "./bwa mem %s %s %s > %s"%(self.refname, self.fq1name, self.fq2name, self.samname)
		samtools_view = "./samtools view -bT %s -o %s %s"%(self.refname, self.bamname, self.samname)
		samtools_sort = "./samtools sort %s %s"%(self.bamname, self.sortbam)
		samtools_index = "./samtools index %s"%(self.sortbam + ".bam")
		subprocess.call(indexing, shell=True)
		subprocess.call(samtools_faidx, shell=True)
		subprocess.call(command, shell=True)
		subprocess.call(samtools_view, shell=True)
		subprocess.call(samtools_sort, shell=True)
		subprocess.call(samtools_index, shell=True)
		
	def runPicard(self):
		picard_command = "java -jar picard.jar AddOrReplaceReadGroups I=%s O=%s RGID=%s RGLB=%s RGPL=ILLUMINA RGPU=%s RGSM=%s"%(self.sortbam + ".bam", self.sortbam_rgadd, self.sample_name, self.sample_name, self.sample_name, self.sample_name)
		subprocess.call(picard_command, shell=True)
		samtools_index = "./samtools index %s"%(self.sortbam_rgadd)
		subprocess.call(samtools_index, shell=True)
		dict_name = self.refname.replace("fa" , "dict")
		if not os.path.isfile(dict_name):
			picard_command2 = "java -jar picard.jar CreateSequenceDictionary REFERENCE=%s OUTPUT=%s"%(self.refname, dict_name)
			subprocess.call(picard_command2, shell=True)
	def run(self):
		self.runMapping()
		self.runPicard()
		
class GATKRunningEngine:
	def __init__(self, refname, bamname, vcfname):
		self.refname = refname
		self.bamname = bamname
		self.vcfname = vcfname
	def runGATK(self):
		command = "java -jar GenomeAnalysisTK.jar -R %s -T HaplotypeCaller -I %s -stand_call_conf 20 -o %s"%(self.refname, self.bamname, self.vcfname)
		subprocess.call(command, shell=True)
	def run(self):
		self.runGATK()
