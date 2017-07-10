#!/usr/bin/python

import re
import os
def get_cigar_num(cigar) :
	return re.split("[A-Z=]+",cigar)[:-1]

def get_cigar_alpha(cigar) :
	return re.split("[0-9]+",cigar)[1:]


class ErrorReporter(object):
	def __init__(self, outpath, prefix, ksize, type='dna'):
		self.outpath = outpath
		self.report_file = self.outpath + "/report/analysis.report"
		self.sam_path = self.outpath + "/sam"
		self.vcf_path = self.outpath + "/vcf"
		self.assembly_path = self.outpath + "/assem"
		self.prefix = prefix
		self.type = type
		self.ksize = ksize
		self.diction = {}
	def check_empty(self, contig_fname):
		fsize = os.stat(contig_fname).st_size
		if fsize == 0 :
			return True
		else:
			return False
	def check_no_seed(self, contig_fname):
		contig_file = open(contig_fname)
		if contig_file.readline().rstrip() == "NO SEED CONTIG":
			contig_file.close()
			return True
		else:
			contig_file.close()
			return False

	def check_broken(self, contig_fname):
		contig_file = open(contig_fname)
		flag = True
		for line in contig_file :
			header = line.rstrip()
			if header.startswith(">"):
				if header.startswith(">BROKEN PATH"):
					node = int(header.split("|")[1].split("-")[0])
					#contig_file.close()
			
				else:
			
					flag = False
		return flag

	def analyse(self, gene, sample_index) :
		contig_fname = self.assembly_path + "/" + self.prefix + "_" + gene + "_" + sample_index + "_k" + str(self.ksize) + ".contig"
		contig_fname_final = self.assembly_path + "/" + self.prefix + "_" + gene + "_" + sample_index + "_k" + str(self.ksize) + ".contig.final.fa"
		int_sample_index = int(sample_index)
		if not self.diction.has_key(gene):
			self.diction[gene] = {int_sample_index : {}}

		else:
			self.diction[gene][int_sample_index] = {}
		if self.check_empty(contig_fname):
			self.diction[gene][int_sample_index] = 'UNKNOWN_ASSEMBLY_FAIL'
			return
		if self.check_empty(contig_fname_final):
			self.diction[gene][int_sample_index] = 'NOPROPER'
			return
		if self.check_no_seed(contig_fname) :
			self.diction[gene][int_sample_index] = 'NOSEED'
			return
		tf = self.check_broken(contig_fname)
		if tf :
			self.diction[gene][int_sample_index] = 'BROKEN'
			return

		sam_file_name = self.sam_path + "/" + gene + "_" + str(sample_index) + "_k" + str(self.ksize) + "." + self.type + ".sam"

		### Now parse sam files
		with open(sam_file_name) as sam:
			for line in sam :
				if line.startswith('@') :
					if "LN:" in line:
						self.length = int(line.split("LN:")[-1])
						continue
					else:
						continue
				data = line.rstrip().split("\t")
				cigar = data[5]
				cigar_nums = get_cigar_num(cigar)
				cigar_alphas = get_cigar_alpha(cigar)
				query = data[0]
				if not self.diction[gene][int_sample_index].has_key(query):
					### Handle multiple lines of SAM file
					self.diction[gene][int_sample_index][query] = ''
				
				if len(cigar_nums) >= 2 :
					self.diction[gene][int_sample_index][query] = 'ERROR'
					continue
				else:
					if self.type == 'dna' :
						match_letter = '='
					elif self.type == 'protein' :
						match_letter = 'M'

					if cigar_alphas[0] == match_letter:
						if int(cigar_nums[0]) == self.length:
							self.diction[gene][int_sample_index][query] = 'ErrorFree'
						else:
							self.diction[gene][int_sample_index][query] = 'ERROR'
					else:
						self.diction[gene][int_sample_index][query] = 'ERROR'

	def write_result(self):
		output = open(self.report_file,"aw")
		if self.type == 'dna' :
			output.write('\n[DNA ERROR FREE ANALYSIS]\n')
		elif self.type == 'protein':
			output.write('\n[PROTEIN ERROR FREE ANALYSIS]\n')
		#print self.diction
		for gene in self.diction:
			output.write(gene)
			sample_index = self.diction[gene].keys()
			sample_index.sort()
			for idx in sample_index:
				if type(self.diction[gene][idx]) == str:
					if self.diction[gene][idx] == 'NOSEED' :
						output.write("\tNOSEED")
						continue
					if self.diction[gene][idx].startswith('BROKEN'):
						output.write('\t%s'%(self.diction[gene][idx]))
						continue
					
					if self.diction[gene][idx] == 'UNKNOWN_ASSEMBLY_FAIL':
						output.write('\tUNKNOWN_ASSEMBLY_FAIL(SEE ALIGNMENT RESULT VIA IGV(Integrative genome viewer) )')
						continue
					if self.diction[gene][idx].startswith('NOPROPER'):
						output.write('\t%s'%(self.diction[gene][idx]))
						continue
				free_check = False
				free_set = []
				for q in self.diction[gene][idx]:
					
					if self.diction[gene][idx][q] == 'ERROR':
						continue
					else:
						free_check = True
						free_set.append(q)

				if free_check :
					synopsis = ["FREE("]
					for entry in free_set:
						synopsis.append(entry)
					synopsis.append(")")
					text = ",".join(synopsis)
					output.write("\t" + text)
				else:
					output.write("\tERROR")
			output.write("\n")
		output.write("\n")
		output.close()
