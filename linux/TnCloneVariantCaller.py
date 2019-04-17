#!/usr/bin/python

"""
TnCloneVariantCaller.py

TnClone's variant caller by parsing SAM file

Author : Thomas Sunghoon Heo
"""

import os
from multiprocessing import Process
from SequenceParser import FAParser
from Sam import Sam

TnCloneVCFHeader="""#TnClone-VCF1.0
#This file is TnClone's variant reporting format
#POS is zero-based coordinate(i.e. the starting base location=0)
#SHORT fields tells if assembled contig is shorter than reference file
#REFID\tPOS\tREF\tALT\tSHORT\n"""

def TnCloneVariatCallSingleThread(field):
	### Reading fasta first
	for f in field:
		sam = f[0]
		vcf = f[1]
		ref = f[2]
		parser = FAParser(ref)
		parser.open()
		fasta = {}
		for id, desc, seq in parser.parse():
			fasta[id] = seq
		parser.close()
		VCF = open(vcf , "w")
		VCF.write(TnCloneVCFHeader)
		SamReader = Sam(sam)
		SamReader.open()
		for record in SamReader.parse_record():
			rname = record.rname
			qname = record.qname
			mapped_pos = record.pos - 1 ## Set as zero base
			alphas = record.cigar_alphas()
			numbers = record.cigar_nums()
			zipped = zip(alphas, numbers)
			ref_seq = fasta[rname]
			query_seq = record.seq
			refloc = 0
			altloc = 0
			if rname == '*':
				continue ### Non-mapped
			if len(ref_seq) > len(query_seq):
				short = True
			else:
				short = False
			for idx, z in enumerate(zipped):
				char = z[0]
				num = z[1]
				if char == 'S':
					if idx == 0 :
						seq = seq[num:]
						
					elif idx == len(zipped) - 1:
						seq = seq[:len(seq) - num]
						
					else:
						nums_before_soft = numbers[:idx]
						nums_after_soft = numbers[idx+1:]
						sum_nums_before = sum(nums_before_soft)
						sum_nums_after = sum(nums_after_soft)
						
						front_seq = seq[:sum_nums_before]
						rear_seq = seq[sum_nums_after:]
						
						seq = front_seq + rear_seq
				if (char == 'M') or (char == '='):
					refloc += num
					altloc += num
				elif char == 'I':
					refbase = '-'
					altbase = query_seq[altloc:altloc + num]
					VCF.write(rname + "\t" + str(refloc) + "\t" + refbase + "\t" + altbase + "\t")
					if short == True:
						VCF.write("YES\n")
					else:
						VCF.write("NO\n")
					altloc += num
				elif char == 'D':
					refbase = ref_seq[refloc : refloc + num]
					altbase = '-'
					VCF.write(rname + "\t" + str(refloc) + "\t" + refbase + "\t" + altbase + "\t")
					if short == True:
						VCF.write("Y\n")
					else:
						VCF.write("N\n")
					refloc += num
				elif char =='X':
					refbase = ref_seq[refloc]
					altbase = query_seq[altloc]
					VCF.write(rname + "\t" + str(refloc) + "\t" + refbase + "\t" + altbase + "\t")
					if short == True:
						VCF.write("Y\n")
					else:
						VCF.write("N\n")
					refloc += num
					altloc += num
				else:
					refloc += num
					altloc += num	
		SamReader.close()
		VCF.close()
def TnCloneVariantCallerParallemMain(sample_info,
					refname,
					k,
					SAM_IN_DIR,
					VCF_OUT_DIR,
					nprocs,
					t="dna"):
	procs = []
	fields = []
	sam_dna_suffix = "_k%d.dna.sam"%(k)
	sam_prot_suffix = "_k%d.protein.sam"%(k)
	dna_suffix = "_k%d.dna.vcf"%(k)
	prot_suffix = "_k%d.protein.vcf"%(k)
	with open(sample_info) as IN:
		for line in IN:
			data = line.rstrip().split()
			sample = data[0]
			N = int(data[1])
			for i in range(1,N+1):
				field = []
				if t == 'dna':
					sam = SAM_IN_DIR + os.sep + sample + "_" + str(i) + sam_dna_suffix
					if os.path.isfile(sam) == False:
						continue
					vcf = VCF_OUT_DIR + os.sep + sample + "_" + str(i) + dna_suffix
				else:
					sam = SAM_IN_DIR + os.sep + sample + "_" + str(i) + sam_prot_suffix
					if os.path.isfile(sam) == False:
						continue
					vcf = VCF_OUT_DIR + os.sep + sample + "_" + str(i) + prot_suffix
				field.append(sam)
				field.append(vcf)
				field.append(refname)
				fields.append(field)
	baseline = len(fields) / nprocs
	remain = len(fields) % nprocs
	for i in range(nprocs):
		if i == nprocs - 1:
			chunk = baseline + remain
			st = baseline * i
			ed = st + chunk
			field = fields[st:ed]
		else:
			chunk = baseline
			st = baseline * i
			ed = st + chunk
			field = fields[st:ed]
		p = Process(target=TnCloneVariatCallSingleThread, args=(field,))
		procs.append(p)

	for p in procs:
		p.start()
	for p in procs:
		p.join()
