#!/usr/bin/python

"""
IndexSorter.py
Sorting with respect to index

Algorithm was first accomplished by Hoon J.
Algorithm was corrected by Hanna S.
All automation process was done by Sunghoon H.

This class code is written by Sunghoon H.

"""
import os
import sys
import gzip
def get_abs_path(path) :
	return os.path.abspath(os.path.dirname(path))

class Sorter(object):

	def __init__(self, input_dir, output_dir, input_config) :
		self.idir = input_dir
		self.odir = output_dir
		self.config = input_config
		self.tn5_specific1 = "AGATGTGTATAAGAGACAG"
		self.tn5_specific2 = "AGATGTGTATAAGAGACAG"
		#self.editor = editor
	def FileOpener(self, fname):
		f = open(fname, "rb")
		if (f.read(2) == "\x1f\x8b") :
			f.seek(0)
			return gzip.GzipFile(fileobj=f)
		else:
			f.seek(0)
			return f
	def sort(self) :
		input_dir = self.idir
		output_dir = self.odir
		tn5_specific1 = self.tn5_specific1
		tn5_specific2 = self.tn5_specific2
		with open(self.config) as config_file :
			for line in config_file :
				data = line.rstrip().split("\t")
				try:
					assert len(data) == 5
				except AssertionError :
					raise ValueError("Error in sorting configuration file. Please check configuration file")

				gene = data[0]
				n_samples = int(data[1])
				tn5_index = data[2]
				fwd_fqfile_names = data[3].split(",")
				rev_fqfile_names = data[4].split(",")
				fwd_fqfiles_full_paths = [ input_dir + "/" + fname for fname in fwd_fqfile_names ]
				rev_fqfiles_full_paths = [ input_dir + "/" + fname for fname in rev_fqfile_names ]

				assert len(fwd_fqfiles_full_paths) == len(rev_fqfiles_full_paths)
				assert len(fwd_fqfiles_full_paths) == n_samples
	
				fq_filenames = zip(fwd_fqfiles_full_paths , rev_fqfiles_full_paths)

				## Make directory
				if not os.path.isdir(output_dir + "/" + gene) :
					os.mkdir(output_dir + "/" + gene)

				## Sorting using Tn5 barcode
				for i , fqfile_set in enumerate(fq_filenames) :
					sample_index = i + 1
					fwd_fqname = fqfile_set[0]
					rev_fqname = fqfile_set[1]
					#self.editor.append("[TnClone::Sorter] Target 1 (%s) : %s\n"%(gene,fwd_fqname))
					#self.editor.append("[TnClone::Sorter] Target 2 (%s) : %s\n"%(gene,rev_fqname))
					#sys.stderr.write("[Target 1 (%s)] : %s\n"%(gene,fwd_fqname))
					#sys.stderr.write("[Target 2 (%s)] : %s\n"%(gene,rev_fqname))
					
					gz_fwd_file = self.FileOpener(fwd_fqname)
					gz_rev_file = self.FileOpener(rev_fqname)
					
					output_fwd_filename = output_dir + "/" + gene + "/" + gene + "_" + str(sample_index) + "_R1.sorted.fastq"
					output_rev_filename = output_dir + "/" + gene + "/" + gene + "_" + str(sample_index) + "_R2.sorted.fastq"
					output_fwd = open(output_fwd_filename , "aw")
					output_rev = open(output_rev_filename , "aw")
					i = 0
					header1 = ""
					header2 = ""
					seq1 = ""
					seq2 = ""
					qual1 = ""
					qual2 = ""
					no_bcd = 0
					for line1 in gz_fwd_file :
						line2 = gz_rev_file.next()
						if i % 4 == 0 :
							header1 = line1.rstrip()
							header2 = line2.rstrip()
							i += 1
						elif i % 4 == 1 :
							seq1 = line1.rstrip()
							seq2 = line2.rstrip()
							i += 1
						elif i % 4 == 2 :
							i += 1 # + synbol
						elif i % 4 == 3 :
							qual1 = line1.rstrip()
							qual2 = line2.rstrip()
	
							### Perform sorting here
							seq1_recover = seq1
							seq2_recover = seq2
							qual1_recover = qual1
							qual2_recover = qual2
							if tn5_specific1 in seq1 and tn5_specific2 in seq2 :
								possible_bcd = seq2.split(tn5_specific1)[0]
								if possible_bcd == tn5_index :
									fwd_loc = seq1.find(tn5_specific1) + len(tn5_specific1)
									rev_loc = seq2.find(tn5_specific1) + len(tn5_specific1)
									seq1_recover = seq1[fwd_loc:]
									seq2_recover = seq2[rev_loc:]
									qual1_recover = qual1[fwd_loc:]
									qual2_recover = qual2[rev_loc:]
									output_fwd.write(header1 + "\n" + seq1_recover + "\n+\n" + qual1_recover + "\n")
									output_rev.write(header2 + "\n" + seq2_recover + "\n+\n" + qual2_recover + "\n")
									i += 1
									# print i
								else:
									no_bcd += 1
									i += 1
									# print i
									continue
							else:
								i += 1
								# print i
								continue
					output_fwd.close()
					output_rev.close()
				
