#!/usr/bin/python

"""
KmerCountingThread.py

Counts k-mers from multiple fastq files in parallel

Required : kmer_counter (executable in C++)

Author : Thomas Sunghoon Heo
"""
import sys
import os
from multiprocessing import Process
import subprocess
cwd = os.getcwd()
def Action(out_path, ksz, fastq_file_fwds, fastq_file_revs, proc_id):
	Z = zip(fastq_file_fwds, fastq_file_revs)
	k = str(ksz)
	sys.stdout.write("[TnClone::KmerCountMultiThread] Process %d in action\n"%(proc_id))
	for z in Z :
		gene = z[0].split("_")[0]
		value = z[0].split("_")[1]
		outfile_name = out_path + os.sep + gene + "_" + value + ".fasta"
		fwd = z[0]
		rev = z[1]
		script = cwd + os.sep + "./kmer_counter %s %s %s %s"%(k, fwd, rev, outfile_name)
		subprocess.call(script,shell=True)




def KmerCountMultiThread(k, trim_path, fastq_1s, fastq_2s, nproc):
	q = len(fastq_1s) / nproc
	r = len(fastq_1s) % nproc
	Ps = []
	for i in range(nproc):
		st = i
		if i == nproc -1 :
			ed = q + r
		else:
			ed = q
		fwds = fastq_1s[st:ed]
		revs = fastq_ws[st:ed]
		P = Process(target=Action, args=(trim_path, k, fwds, revs,))
		Ps.append(P)

	for p in Ps:
		p.start()

	for p in Ps:
		P.join()
