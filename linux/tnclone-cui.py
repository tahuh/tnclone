#!/usr/bin/python

"""
tnclone-cui.py
The command line interface version of TnClone
Use this code when you are expert at Linux system!
Many of classes are imported from GUI version of TnClone
This source can be used at any Linux system with python installed

This code will be used for future Web based TnClone
Author : Sunghoon Heo
"""

### Import libraries and classes
### Python built-in libraries
import sys
import os
import subprocess
import argparse
import logging
import logging.handlers
import time

### Class imports
from IndexSorter import Sorter
from TnCloneSorter import TnCloneSortParallelMain
from Trimmer import Trimmer
from Sampler import Sampler
from SequenceParser import FQParser, FAParser
from Graph import Graph
from KmerCounter import KmerCounter
#from Modify import MinimumDepthKmerCollector
from Modifier import Modifier2
from Assembly import AssemblyEngine
from AlignmentManager import SSWAlignmentManager
from Analytics import ReportAnalysis , ConfidentSeqSelector
### Argument parser
parser = argparse.ArgumentParser(description="TnClone - Command Line User Interface")

### Required fields
parser.add_argument("--sample_file" , type=str, required=True, help="Sample information file")
parser.add_argument("--output_dir", type=str, required=True, help="Output directory")

### Non-required fields
parser.add_argument("--ngs_directory", type=str, default=None, help="Raw NGS directory")
parser.add_argument("--sort_info", type=str, default=None, help="Information for sample sorting")
parser.add_argument("--start_seq", type=str, default=None, help="Start of assembly sequence")
parser.add_argument("--end_seq", type=str, default=None, help="End of assembly sequence")
parser.add_argument("--reference_folder", type=str, default=None, help="A directory that contains set of reference files")
parser.add_argument("--reference_file", type=str, default=None, help="Reference FILE which contains multiple references")
parser.add_argument("--region_file", type=str, default=None, help="Assembly region information file. Must use without --start_seq/--end_seq")
parser.add_argument("--kmer_size", type=int, default=63, help="Length of k-mer. default=63")
parser.add_argument("--min_kmer_occurence", type=int, default=3, help="Minimum k-mer occurence. default=3")
parser.add_argument("--downsample_rate", type=float, default=0.3, help="Rate of downsampling. default=0.3")
parser.add_argument("--number_of_seed_mismatch_allowed", type=int, default=3, help="Number of mismatch bases allowed for alterlative seed search. default=1")
parser.add_argument("--alignment_match_score", type=int, default=1, help="Alignment match score for SSW. default=1")
parser.add_argument("--alignment_mismatch_penalty", type=int, default=3, help="Alignment mismatch penalty. must be positive. default=3")
parser.add_argument("--alignmemt_gap_open_penalty", type=int, default=5, help="Alignment gap open penalty. must by positive. default=5")
parser.add_argument("--alignment_gap_ext_penalty", type=int, default=3, help="Alignment gap extension penalty. must be positive. default=3")
parser.add_argument("--node_to_node_ratio", type=float, default=0.01, help="Minimum k-mer ratio between two nodes of an edge. default=0.01")
parser.add_argument("--multiple_edge_ratio", type=float, default=0.05, help="Minimum k-mer ratio between multiple edges. default=0.05")
parser.add_argument("--alternative_seed_selection_method", type=str, default="mismatch", help="Alternative seed selection method(mismatch/depth). default=mismatch")
parser.add_argument("--num_cores", type=int, default=4 ,help="Number of physical cores to use. default=4")
### Argument controller actions
parser.add_argument("--denovo_assembly", default=False, action="store_true", help="Perform de novo assembly only. Reference file/folder MUST NOT set")
parser.add_argument("--no_confident_seq", default=False, action="store_true", help="Not extracting confident sequences if set. default=OFF")
parser.add_argument("--perform_downsample", default=False, action="store_true", help="Whether to perform downsampling defore assembly to speed up")
parser.add_argument("--use_alternative_seed", default=False, action="store_true", help="Force select alternative seed sequence")
parser.add_argument("--no_sort", default=False, action="store_true", help="Not perform sample sorting")
parser.add_argument("--no_trim", default=False, action="store_true", help="Not perform trimming")
parser.add_argument("--no_assembly", default=False,action="store_true", help="Not perform assembly")
parser.add_argument("--no_analysis", default=False, action="store_true", help="Not perform downstream analysis")

TNCLONE_ARGS = parser.parse_args()

### Setup logger
logger = logging.getLogger("TnCloneLogger")
formatter = logging.Formatter("[%(levelname)s|%(filename)s:%(lineno)s] %(asctime)s > %(message)s")
fileHandler = logging.FileHandler("TnCloneLogFile-" + str(time.time()) + ".log")
streamHandler = logging.StreamHandler()
fileHandler.setFormatter(formatter)
streamHandler.setFormatter(formatter)
logger.addHandler(fileHandler)
logger.addHandler(streamHandler)
logger.setLevel(logging.DEBUG)

tnclone_logo = """ _________           ______  __                         
|  _   _  |        .' ___  |[  |                        
|_/ | | \_|_ .--. / .'   \_| | |  .--.   _ .--.  .---.  
    | |   [ `.-. || |        | |/ .'`\ \[ `.-. |/ /__\\ 
   _| |_   | | | |\ `.___.'\ | || \__. | | | | || \__., 
  |_____| [___||__]`.____ .'[___]'.__.' [___||__]'.__.' """

if sys.version_info[0] == 3:
	print(tnclone_logo)
else:
	print tnclone_logo
	
def PrintLog(msg):
	logger.info(msg)

def WarnLog(msg):
	logger.warning(msg)

def CriticalLog(msg):
	logger.critical(msg)

def ErrorLog(msg):
	logger.error(msg)

def ShowParameters():
	PrintLog("TnClone argument is shown below")
	s = 'tnclone-cui.py'
	for arg in vars(TNCLONE_ARGS):
		attr = getattr(TNCLONE_ARGS, arg)
		tmp = " --" + str(arg) + "=" + str(attr)
		s += tmp
	PrintLog(s)
def LaunchCommand(command):
	subprocess.call(command, shell=True)
def ValidateInputArguments():
	global TNCLONE_ARGS
	args = TNCLONE_ARGS
	if args.no_sort == False:
		if args.sort_info == None:
			ErrorLog("--no_sort is set False but no --sort_info. Please check")
			exit(-1)
		else:
			if args.ngs_directory == None:
				ErrorLog("--no_sort is set False but no NGS data. Please check")
				exit(-1)
	if args.no_assembly == False:
		if args.start_seq != None and args.end_seq == None:
			ErrorLog("--start_seq is set but not --end_seq. Please set --end_seq")
			exit(-1)
		if args.start_seq == None and args.end_seq != None:
			ErrorLog("--end_seq is set but not --start-seq. Please set --start_seq for assembly")
			exit(-1)
		if args.start_seq == None and args.end_seq == None:
			if args.region_file == None:
				ErrorLog("--start_seq and --end_seq has empty value but also --region_file is empty. Please set")
				exit(-1)
		if args.region_file != None:
			if args.start_seq != None or args.end_seq != None:
				ErrorLog("--region_file is set but either --start_seq or --end_seq if set.")
				exit(-1)
			else:
				if args.refernce_file == None and args.reference_file != None:
					ErrorLog("--region_file is set but no reference information is given")
				exit(-1)
	if args.denovo_assembly:
		if args.reference_folder != None or args.reference_file != None:
			WarnLog("--denovo_asasembly is set but also reference is given. Force --no_analysis")
			TNCLONE_ARGS.no_analysis = True
			TNCLONE_ARGS.reference_file = None
			TNCLONE_ARGS.reference_folder = None

def driver():
	SORT_PATH = TNCLONE_ARGS.output_dir + os.sep + "sorted"
	TRIM_PATH = TNCLONE_ARGS.output_dir + os.sep + "trimmed"
	if TNCLONE_ARGS.no_sort == False:
		PrintLog("De-multiplex samples")
		if not os.path.isdir(TNCLONE_ARGS.output_dir + os.sep + "sorted"):
			os.mkdir(TNCLONE_ARGS.output_dir + os.sep + "sorted")
		#sorter = Sorter(TNCLONE_ARGS.ngs_directory, SORT_PATH, TNCLONE_ARGS.sort_info)
		#sorter.sort()
		TnCloneSortParallelMain(TNCLONE_ARGS.sort_info, "AGATGTGTATAAGAGACAG", TNCLONE_ARGS.output_dir + os.sep + "sorted", 4)
		PrintLog("Done De-multiplexing")
	if TNCLONE_ARGS.no_trim == False:
		if not os.path.isdir(TNCLONE_ARGS.output_dir + os.sep + "trimmed"):
			os.mkdir(TNCLONE_ARGS.output_dir + os.sep + "trimmed")
		PrintLog("Trimming using Trimmomatic")
		trimmer = Trimmer(SORT_PATH, TRIM_PATH, TNCLONE_ARGS.sample_file, is_cmd = True)
		trimmer.trim()
	if TNCLONE_ARGS.no_assembly == False:
		### Count-up k-mers
		INIT_PATH = TRIM_PATH
		with open(TNCLONE_ARGS.sample_file) as FILE:
			for line in FILE:
				data = line.rstrip().split()
				n_samples = int(data[1])
				sample_name = data[0]
				for i in range(1,n_samples+1):
					PREFIX = INIT_PATH + os.sep + sample_name + os.sep + sample_name + "_" + str(i) + "_R"
					F1 = PREFIX + "1.trimmed.fastq"
					F2 = PREFIX + "2.trimmed.fastq"
					if TNCLONE_ARGS.perform_downsample:
						sampler = Sampler(F1,F2, sample_ratio = TNCLONE_ARGS.downsample_rate)
						sampler.sample()
						PrintLog("Down sample is set with rate %f. Perform down sample."%(TNCLONE_ARGS.downsample_rate))
						F1 = PREFIX + "1.trimmed.fastq.dwn.fastq"
						F2 = PREFIX + "2.trimmed.fastq.dwn.fastq"
					
					P1 = FQParser(F1)
					P2 = FQParser(F2)
					P1.open(); P2.open()
					
					#kmer_counter = KmerCounter(P1, P2, TNCLONE_ARGS.kmer_size)
					#kmer_counter.Count()
					
					#PrintLog("Selecting k-mers lower than minimum occurence set as %d"%(TNCLONE_ARGS.min_kmer_occurence))
					#unused_kmer_collector = MinimumDepthKmerCollector(kmer_counter, TNCLONE_ARGS.min_kmer_occurence)
					#unused_kmer_collector.Action()
					PrintLog("Execute for %d-th sample"%(i))
					PrintLog("Build naive de Bruijn graph before assembly")
					
					P1.rewind()
					P2.rewind()
					
					graph = Graph(P1,P2, TNCLONE_ARGS.kmer_size)
					
					graph.build_naive()
					
					PrintLog("Done building naive de Bruijn graph")
					
					PrintLog("Modifying naive de Bruijn graph")
					initial_num_nodes = graph.num_nodes()
					PrintLog("%d nodes initially present"%(initial_num_nodes))
					# modification_engine = GraphModificationEngine(graph, unused_kmer_collector, TNCLONE_ARGS.node_to_node_ratio, TNCLONE_ARGS.multiple_edge_ratio)
					# modification_engine.Action()
					
					# graph= modification_engine.ShowGraph()
					modifier = Modifier2(graph, TNCLONE_ARGS.min_kmer_occurence, edge_freq = TNCLONE_ARGS.node_to_node_ratio, edge_freq2 = TNCLONE_ARGS.multiple_edge_ratio)
					try:
						modifier.enforce_modify()
					except ValueError:
						WarnLog("TnClone cannot undergo assembly for sample %s due to errorneous k-mers"%(data[0] + "_"+ str(i)))
						WarnLog("Island: %d, Dead end: %d, Low Depth: %d, Low Freq: %d"%(modifier.island, modifier.dead_end, modifier.low_depth, modifier.low_allele))
						continue
					
					graph = modifier.graph
					#del unused_kmer_collector ### Explicitly tell GC to work!
					
					modified_num_nodes = graph.num_nodes()
					
					PrintLog("Done modification(%d nodes -> %d nodes). We got ready to assemble graph"%(initial_num_nodes, modified_num_nodes))
					
					if modified_num_nodes == 0 :
						WarnLog("TnClone cannot assemble since there is no nodes")
						tmp1 = TNCLONE_ARGS.output_dir + os.sep + "assem" + os.sep + sample_name + "_" + str(i) + "_k" + str(TNCLONE_ARGS.kmer_size) + ".contig"
						with open(tmp1 , "w") as F: pass
						tmp1 = TNCLONE_ARGS.output_dir + os.sep + "assem" + os.sep + sample_name + "_" + str(i) + "_k" + str(TNCLONE_ARGS.kmer_size) + ".contig.indel"
						with open(tmp1 , "w") as F: pass
						tmp1 = TNCLONE_ARGS.output_dir + os.sep + "assem" + os.sep + sample_name + "_" + str(i) + "_k" + str(TNCLONE_ARGS.kmer_size) + ".contig.gool"
						with open(tmp1 , "w") as F: pass
						tmp1 = TNCLONE_ARGS.output_dir + os.sep + "assem" + os.sep + sample_name + "_" + str(i) + "_k" + str(TNCLONE_ARGS.kmer_size) + ".contig.final.fa"
						with open(tmp1 , "w") as F: pass
						continue
					PrintLog("Start assembly")
					PrintLog("Seed: %s, Terminal: %s"%(TNCLONE_ARGS.start_seq.upper(), TNCLONE_ARGS.end_seq.upper()))
					assembler = AssemblyEngine(TNCLONE_ARGS.start_seq.upper() , TNCLONE_ARGS.end_seq.upper(), graph)
					paths = assembler.Action()
					paths = list(paths)
					
					FILE_PATH = TNCLONE_ARGS.output_dir + os.sep + "assem" + os.sep
					if not os.path.isdir(FILE_PATH):
						os.mkdir(FILE_PATH)
					FILE_PATH = FILE_PATH + sample_name + "_" + str(i) + "_k" + str(TNCLONE_ARGS.kmer_size) + ".contig"
					STREAM = open(FILE_PATH, "w")
					tmpi = 0
					for p in paths:
						seq = p[0]
						for pp in p[1:]:
							seq += pp[-1]
						tmpi+=1
						STREAM.write(">" + str(tmpi) + "\n" + seq + "\n")
					STREAM.close()
					if os.stat(FILE_PATH).st_size == 0 :
						tmp1 = TNCLONE_ARGS.output_dir + os.sep + "assem" + os.sep + sample_name + "_" + str(i) + "_k" + str(TNCLONE_ARGS.kmer_size) + ".contig.indel"
						with open(tmp1 , "w") as F: pass
						tmp1 = TNCLONE_ARGS.output_dir + os.sep + "assem" + os.sep + sample_name + "_" + str(i) + "_k" + str(TNCLONE_ARGS.kmer_size) + ".contig.gool"
						with open(tmp1 , "w") as F: pass
						tmp1 = TNCLONE_ARGS.output_dir + os.sep + "assem" + os.sep + sample_name + "_" + str(i) + "_k" + str(TNCLONE_ARGS.kmer_size) + ".contig.final.fa"
						with open(tmp1 , "w") as F: pass
	if TNCLONE_ARGS.no_confident_seq == False:
		#if TNCLONE_ARGS.denovo_assembly == True:
		with open(TNCLONE_ARGS.sample_file) as F:
			for line in F:
				data = line.rstrip().split("\t")
				sample_name = data[0]
				num = int(data[1])
				for i in range(1, num+1):
					FILE_PATH = TNCLONE_ARGS.output_dir + os.sep + "assem" + os.sep
					infile_name = FILE_PATH + sample_name + "_" + str(i) + "_k" + str(TNCLONE_ARGS.kmer_size) + ".contig"
					try:
						infile = open(infile_name)
					except IOError:
						tmp1 = TNCLONE_ARGS.output_dir + os.sep + "assem" + os.sep + sample_name + "_" + str(i) + "_k" + str(TNCLONE_ARGS.kmer_size) + ".contig.final.fa"
						with open(tmp1 , "w") as F: 
							pass
						continue
					if os.stat(infile_name).st_size == 0 :
						tmp1 = TNCLONE_ARGS.output_dir + os.sep + "assem" + os.sep + sample_name + "_" + str(i) + "_k" + str(TNCLONE_ARGS.kmer_size) + ".contig.final.fa"
						with open(tmp1 , "w") as F: 
							pass
						continue
					cwd = os.getcwd()
					MAFFT_PATH = cwd + os.sep + "mafft"
					MAFFT_OUT_PATH = TNCLONE_ARGS.output_dir + os.sep + "sam"
					if not os.path.isdir(MAFFT_OUT_PATH):
						os.mkdir(MAFFT_OUT_PATH)
					MAFFT_OUT = TNCLONE_ARGS.output_dir + os.sep + "sam" + os.sep + sample_name + "_k" + str(TNCLONE_ARGS.kmer_size) + ".mafft"
					MAFFT_COMMAND = MAFFT_PATH + " " + infile_name + " > " + MAFFT_OUT
					LaunchCommand(MAFFT_COMMAND)
					
					if not os.path.isdir(TNCLONE_ARGS.output_dir + os.sep + "vcf"):
						os.mkdir(TNCLONE_ARGS.output_dir + os.sep + "vcf")
					MSA2VCF_OUTPUT = TNCLONE_ARGS.output_dir + os.sep + "vcf" + os.sep + sample_name + "_k" + str(TNCLONE_ARGS.kmer_size) + ".dna.vcf"
					MSA2VCF_COMMAND = "java -jar " + cwd + os.sep + "msa2vcf.jar " + MAFFT_OUT + " > " + MSA2VCF_OUTPUT
					LaunchCommand(MSA2VCF_COMMAND)
					
					
					PREFIX = TRIM_PATH + os.sep + sample_name + os.sep + sample_name + "_" + str(i) + "_R"
					F1 = PREFIX + "1.trimmed.fastq"
					F2 = PREFIX + "2.trimmed.fastq"
					
					
					CSS = ConfidentSeqSelector(infile_name,MSA2VCF_OUTPUT,F1,F2,TNCLONE_ARGS.kmer_size)
					try:
						CSS.Open()
					except IOError:
						continue
					PrintLog("Build k-mer depth table for confident contig selection")
					CSS.BuildDepthTable()
					PrintLog("Read FASTA file")
					CSS.ReadFASTA()
					PrintLog("Reading VCF")
					CSS.ReadVCFInfo()
					PrintLog("Extract confident sequence")
					CSS.Extract()
					PrintLog("Flush output")
					tmp_fname = TNCLONE_ARGS.output_dir + os.sep + "assem" + os.sep + sample_name + "_" + str(i) + "_k" + str(TNCLONE_ARGS.kmer_size) + ".contig.final.fa"
					CSS.FlushSelected(tmp_fname)
					CSS.Close()
	if TNCLONE_ARGS.no_analysis == False:
					PrintLog("Confident assembled contig selection step initiated")
					if not os.path.isdir(TNCLONE_ARGS.output_dir + os.sep + "pre-aln"):
						PrintLog("Generate pre-aln directory to select non-indel DNA sequences assembled\n")
						os.mkdir(TNCLONE_ARGS.output_dir + os.sep + "pre-aln")
						
					PrintLog("Analyse assembled contigs whether those have indels or not")
					if TNCLONE_ARGS.reference_file != None:
						aligner = SSWAlignmentManager( TNCLONE_ARGS.alignment_match_score,
														TNCLONE_ARGS.alignment_mismatch_penalty,
														TNCLONE_ARGS.alignmemt_gap_open_penalty,
														TNCLONE_ARGS.alignment_gap_ext_penalty,
														ref_file = TNCLONE_ARGS.reference_file,
														ref_dir = None)
						PRE_ALN_INPUT_FILE = FILE_PATH 
						PRE_ALN_OUTPUT_FILE = TNCLONE_ARGS.output_dir + os.sep + "pre-aln" + os.sep + sample_name + "_" + str(i)  + "_k" + str(TNCLONE_ARGS.kmer_size) + "_pre-aln.dna.sam"
						aligner.PreAlignActionDNA(PRE_ALN_OUTPUT_FILE,PRE_ALN_INPUT_FILE,REFERENCE=None)
					else:
						aligner = SSWAlignmentManager( TNCLONE_ARGS.alignment_match_score,
														TNCLONE_ARGS.alignment_mismatch_penalty,
														TNCLONE_ARGS.alignmemt_gap_open_penalty,
														TNCLONE_ARGS.alignment_gap_ext_penalty,
														ref_file = None,
														ref_dir = TNCLONE_ARGS.reference_folder)
						PRE_ALN_INPUT_FILE = FILE_PATH 
						PRE_ALN_OUTPUT_FILE = TNCLONE_ARGS.output_dir + os.sep + "pre-aln" + os.sep + sample_name + "_" + str(i)  + "_k" + str(TNCLONE_ARGS.kmer_size) + "_pre-aln.dna.sam"
						REF_NAME = TNCLONE_ARGS.reference_folder + os.sep + sample_name + ".gene.fa"
						aligner.PreAlignActionDNA(PRE_ALN_OUTPUT_FILE,PRE_ALN_INPUT_FILE,REFERENCE=REF_NAME)
													
					PrintLog("Done alignment")
if __name__ == "__main__":
	PrintLog("Validate TnClone input parameters")
	ValidateInputArguments()
	ShowParameters()
	driver()
