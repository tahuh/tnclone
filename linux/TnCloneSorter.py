#!/usr/bin/python

"""
TnCloneSorter.py
TnClone now use parallel file sorting
Each thread handles one line of config file each

Author : Sunghoon Heo
"""
import os
import gzip
from multiprocessing import Process

def format_checker(fqname) :
	F = open(fqname,"rb")
	bytes = F.read(3)
	if bytes == "\x1f\x8b\x08" :
		F.seek(0)
		return gzip.GzipFile(fqname , fileobj=F)
	else:
		F.seek(0)
		return F

class TnCloneSorterSingleProcessingUnit:
	def __init__(self,infile1, infile2, outfile1, outfile2, me, fwd, rev):
		self.infile1 = format_checker(infile1)
		self.infile2 = format_checker(infile2)
		self.outfile1 = open(outfile1, "aw")
		self.outfile2 = open(outfile2, "aw")
		self.me = me
		self.fwd = fwd
		self.rev = rev

	def Action(self):
		i = 0
		id1 = ''
		id2 = ''
		seq1 = ''
		seq2 = ''
		qual1 = ''
		qual2 = ''
		for line1 in self.infile1:
			line2 = self.infile2.next()
			if i % 4 == 0:
				id1 = line1.rstrip()
				id2 = line2.rstrip()
				i += 1
				continue
			elif i % 4 == 1:
				seq1 = line1.rstrip()
				seq2 = line2.rstrip()
				i += 1
				continue
			elif i % 4 == 2:
				i += 1
				continue
			elif i % 4 == 3:
				qual1 = line1.rstrip()
				qual2 = line2.rstrip()

				### Search ME sequence
				if (self.me not in seq1) or (self.me not in seq2):
					continue

				### Search for tn5 barcodes
				fwdi = seq1.find(self.fwd)
				revi = seq2.find(self.rev)

				if fwdi == -1 or revi == -1:
					continue
				
				fwd_bcd = seq1.split(self.me)[0]
				rev_bcd = seq2.split(self.me)[0]
				if (fwd_bcd == self.fwd) and ( rev_bcd == self.rev ):
					fwd_me_loc = self.seq1.find(self.me)
					rev_me_loc = self.seq2.find(self.me)
					seq1 = seq1[fwd_me_loc + len(self.me):]
					seq2 = seq2[fwd_me_loc + len(self.me):]
					qual1 = seq1[fwd_me_loc + len(self.me):]
					qual2 = seq2[fwd_me_loc + len(self.me):]
					self.outfile1.write(id1 + "\n" + seq1 + "\n+\n" + qual1 + "\n")
					self.outfile2.write(id2 + "\n" + seq2 + "\n+\n" + qual2 + "\n")

		self.outfile1.close()
		self.outfile2.close()

def TnCloneSortFunction(data, me, odir):
	for line in data:
		splits = line.rstrip().split()
		gene = splits[0]
		num = int(splits[1])
		fwd = splits[2]
		rev = splits[3]
		ngs_list1 = splits[4].split(",")
		ngs_list2 = splits[5].split(",")

		for i in range(1, num+1):
			infile1 = ngs_list1[i-1]
			infile2 = ngs_list2[i-1]
			try:
				os.mkdir(odir + os.sep + gene)
			except OSError:
				pass
			outfile1 = odir + os.sep + gene + os.sep + gene + "_" + str(i) + "_R1.sorted.fastq"
			outfile2 = odir + os.sep + gene + os.sep + gene + "_" + str(i) + "_R2.sorted.fastq"
			SortEngine = TnCloneSorterSingleProcessingUnit(infile1, infile2, outfile1, outfile2, me, fwd, rev)
			SortEngine.Action()

def TnCloneSortParallelMain(config_file, me, output_dir, nt):
	F = open(config_file)
	file_lines = list(F.readlines())
	N = len(file_lines)
	chunk = int(N / nt)
	remain = N % nt
	num_data = 0
	procs = []
	for i in range(nt):
		if i == nt - 1:
			num_data = chunk + remain
		else:
			num_data = chunk

		data_pass = file_lines[i * chunk : i * chunk + num_data]
		p = Process(target=TnCloneSortFunction, args=(data_pass, me, output_dir))
		procs.append(p)

	for i in range(nt):
		procs[i].start()

	for i in range(nt):
		procs[i].join()
