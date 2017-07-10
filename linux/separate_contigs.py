#!/usr/bin/python

from argparse import ArgumentParser

try:
	from Bio import SeqIO as SeqParser
	BIOPY_AVAIL = True
except ImportError:
	from SequenceParser import FAParser as SeqParser
	BIOPY_AVAIL = False

parser = ArgumentParser()
parser.add_argument("--contig_file" , type=str, required=True)
parser.add_argument("--node_file" , type=str,required=True)
parser.add_argument("--path_score" , type=float,required=False)
parser.add_argument("--cv_limit" , type=float, required=False)
parser.add_argument("--ksize" , type=int,required=True)

args = parser.parse_args()
contig_file = open(args.contig_file)
node_file = open(args.node_file)
ksize = args.ksize

path_score_limit = args.path_score
cv_limit = args.cv_limit


#### BUILD TABLE

index_table = {}
depth_table = {}
for line in node_file :
	data = line.rstrip().split("\t")
	curnode = data[0]
	curnode_depth = int(data[1])
	if not index_table.has_key(curnode):
		index_table[curnode] = True
	if not depth_table.has_key(curnode):
		depth_table[curnode] = curnode_depth
	if len(data) > 2 :
		next_nodes = [node[1:] + b for b in data[2].split(",")]
		next_depths = [int(x) for x in data[3].split(",")]
		for i,n in enumerate(next_nodes) :
			if not index_table.has_key(n) :
				index_table[n] = True
			if not depth_table.has_key(n):
				depth_table[n] = next_depths[i] 


### Indexing
node_set = index_table.keys()
node_set.sort()

index = 0
inverse_table = {}
for n in node_set:
	index_table[n] = index
	inverse_table[index] = n
	index += 1

### Retreiving sequence
seq_dict = {}

if BIOPY_AVAIL :
	for record in SeqIO.parse(cocntig_file, "fasta"):
		ID = record.seq
		seq = str(record.seq)
		seq_dict[ID] = seq

else:
	parser = FAParser(args.contig_file)
	parser.open()
	for id , seq in parser.parse() :
		seq_dict[id] = seq

### Then k-mersize and insert into differenct list using indexes
num_contigs = len(seq_dict.keys())

seq_ids = seq_dict[keys]

index_set = [ [None] for _ in xrange(len(seq_ids)) ]

for i , val in enumerate(seq_ids) :
	seq = seq_dict[val]
	for j in xrange(len(seq)-ksize):
		kmer = seq[j:j+ksize]
		idx = index_table[kmer]
		index_set[i].append(idx)
