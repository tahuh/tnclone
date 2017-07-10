#!/usr/bin/python

"""
This code automatically trims ME sequence only if one used single sample

Author : Sunghoon Heo
"""

import os
from argparse import ArgumentParser

from memory import MemoryUsageManager
from SequenceParser import FQParser

from Bio import SeqIO
import gzip
import sys

pid = os.getpid()

manager = MemoryUsageManager(pid , unit = "Gb")

fwd_fq_parser = FQParser()
rev_fq_parser = FQParser()

MESEQ = "AGATGTGTATAAGAGACAG"

"""
File format
<gene><tab><fwd_fq><tab><rev_fq>
"""

ArgParser = ArgumentParser(description="Single sample Trimmer")

ArgParser.add_argument("--input_config" , type=str,required=True, help="Input configuration")
ArgParser.add_argument("--output_dir" , type=str,required=True,help="Output directory")
ArgParser.add_argument("--raw_path" , type=str, required=True , help="RAW NGS PATH")

args = ArgParser.parse_args()

def abspath(path) :
	return os.path.abspath(path)


input_config_file = args.input_config
outdir = abspath(args.output_dir)
rawpath = abspath(args.raw_path)

### Reading config file

conf_dict = {}
with open(input_config_file) as config_file:
	for line in config_file:
		data = line.rstrip().split("\t")
		gene_id = data[0]
		fwd_fqname = data[1]
		rev_fqname = data[2]
		if not conf_dict.has_key(gene_id) :
			conf_dict[gene_id] = {"fwd" : [fwd_fqname] , "rev" : [rev_fqname] ,"items" : 1}
		else:
			conf_dict[gene_id]["fwd"].append(fwd_fqname)
			conf_dict[gene_id]["rev"].append(rev_fqname)
			conf_dict[gene_id]["items"] += 1

### Process to trim sequence
SAMPLE_OUT_SUFFIX_FWD = "_R1.sorted.fastq"
SAMPLE_OUT_SUFFIX_REV = "_R2.sorted.fastq"

for gene in conf_dict.keys():
	n_genes = conf_dict["items"]
	for x in xrange(n_genes) :
		sample_out_name_prefix = gene + "_" + str(x+1)
		fwd_oname = sample_out_name_prefix + SAMPLE_OUT_SUFFIX_FWD
		rev_oname = sample_out_name_prefix + SAMPLE_OUT_SUFFIX_REV
		
		fwd_out = open(fwd_oname , "aw")
		rev_out = open(rev_oname , "aw")

		fq_fwd_name = conf_dict["fwd"][x]
		fq_rev_name = conf_dict["rev"][x]

		#fwd_fq_parser.open(fq_fwd_name)
		#rev_fq_parser.open(fq_rev_name)

		fwd_fq = gzip.open(fq_fwd_name)
		rev_fq = gzip.open(fq_rev_name)

		header1 = ""
		header2 = ""
		seq1 = ""
		seq2 = ""
		qual1 = ""
		qual2 = ""

		for i , line in enumerate(fwd_fq) :
			line2 = rev_fq.next()

			if i & 3 == 0 :
				header1 = line.rstrip()
				header2 = line2.rstrip()
			elif i & 3 == 1 :
				seq1 = line.rstrip()
				seq2 = line.rstrip()
			elif i & 3 == 2 :
				continue
			elif i & 3 == 3:
				qual1 = line.rstrip()
				qual2 = line.rstrip()

				if MESEQ in seq1 and MESEQ in seq2 :
					pos1 = seq1.find(MESEQ)
					pos2 = seq2.find(MESEQ)
					seq1 = seq1.split(MESEQ)[1]
					seq2 = seq2.split(MESEQ)[1]
					qual1 = qual1[pos1 + len(MESEQ):]
					qual2 = qual2[pos1 + len(MESEQ):]

					fwd_out.write(header1 + "\n" + seq1 + "\n+\n" + qual1 + "\n")
					rev_out.write(header2 + "\n" + seq2 + "\n+\n" + qual2 + "\n")

		fwd_out.close()
		rev_out.close()

manager.profile()

sys.stderr.write("[MEMPROF] %s"%(manager.mem_usage))
