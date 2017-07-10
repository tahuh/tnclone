#!/usr/bin/python

"""
tnclone_cmd.py

Linux CUI version of TnClone

This script will not create GUI for your analysis

Author : Sunghoon Heo
"""



import time
import os
import sys
import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import optparse

src_path = os.path.dirname(os.path.abspath(__file__))


### CUSTOM module (Algorithms)
from IndexSorter import Sorter
from Trimmer import Trimmer
from KFS import KFS
from FastQParser import SequencingReadParser as SRP
from AssemblyMachine2 import Assembler
from SrcSinkFinder import SrcSinkFinder as SSF
from DupRemoval import DupRemovalTool as DRT
from FinalReporter import AminoAcidReporter , DNAReporter, ContigNumberReporter
from memory import MemoryUsageManager
from Stats import ContigStats2
from ErrorSeedCorrector import ErrorSeedCorrector
from Sampler import Sampler
from ReportError import ErrorReporter
try:
	from Bio import SeqIO
	from Bio.Seq import Seq
	seqio_avail = True
except ImportError:
	from CustomParser import FastaParser
	from DNAUtils import translate
	seqio_avail = False
	
AUTHOR = "Sunghoon Heo"
AUTHOR_EMAIL = "team.tnclone@gmail.com"

HELP_MESSAGE="""Usage: tnclone-cmd.py [options]

TnClone CUI version

Options:
  -h, --help            show this help message and exit
  --ngs_path=NGS_PATH   Raw NGS result path
  --sort_file=SORT_FILE
                        Sorting information containing file.
  --sample_file=SAMPLE_FILE
                        Sample information file.
  --start=START         Start sequence of assembly. Length must be same as
                        k-mer size option
  --end=END             End sequence of assembly. Length must be same as k-mer
                        size option
  --ref_path=REF_PATH   Reference directory that contains reference files for
                        analysis
  --output_path=OUTPUT_PATH
                        Output path
  --kmer_size=KMER_SIZE
                        K-mer size. default = 63
  --min_kocc=MIN_KOCC   Minimum k-mer occurence.
                        This value is used to reduce ambiguity. default = 3
  --sort_on=SORT_ON     Check whether to turn on sorting. 0 for no 1 for yes.
                        default = 1
  --trim_on=TRIM_ON     Check whether to turn on trimming. 0 for no 1 for yes.
                        default = 1
  --assembly_on=ASSEMBLY_ON
                        Check whether to turn on assembly. 0 for no 1 for yes.
                        default = 1
  --analysis_on=ANALYSIS_ON
                        Check whether to turn on analysis. 0 for no 1 for yes.
                        default = 1
  --auto_analysis=AUTO_ANALYSIS
                        Generates analysis automatically. Will start from
                        analysis you specified
                        0 for no 1 for yes. default = 1
  --bed_file=BED_FILE   Bed file
  --down_sampling=DOWN_SAMPLING
                        Downsampling or not 1 for yes and 0 for not. 
                        default = 0
  --down_sample_rate=DOWN_SAMPLE_RATE
                        Down sampling ratio. default = 0.3
  --seed_mismatch_correction_method=SEED_MISMATCH_CORRECTION_METHOD
                        Type of methode to correct seed error. One of depth or
                        mismatch.
                        default=depth
  --num_mismatch=NUM_MISMATCH
                        Number of mismatch allowed.
                        if seed mismatch correction method is set depth , this
                        will not invoked.
                        default = 3
  --gpu_avail=GPU_AVAIL
                        Let the program know if gpu (NVidia card) avail.
                        0 for no 1 for yes. default = 0
  --aln_match=ALN_MATCH
                        Match score. INT value. default = 1
  --aln_mismatch=ALN_MISMATCH
                        Mismatch score. INT value. Will work as negative
                        value. default = 3
  --aln_open_penalty=ALN_OPEN_PENALTY
                        Gap open penalty. INT value. default = 5
  --aln_gap_extension=ALN_GAP_EXTENSION
                        Gap extension penalty. INT value. Default = 3
  --draw_coverage=DRAW_COVERAGE
                        Draw coverage plot before analyis. 0 for no 1 for yes.
                        default = 0"""
						
def show_help():
	print HELP_MESSAGE
	
def show_header_program():

	tnclone_logo = """___________    _________ .__                        
\__    ___/___ \_   ___ \|  |   ____   ____   ____  
  |    | /    \/    \  \/|  |  /  _ \ /    \_/ __ \ 
  |    ||   |  \     \___|  |_(  <_> )   |  \  ___/ 
  |____||___|  /\______  /____/\____/|___|  /\___  >
             \/        \/                 \/     \/ """
	
	print "///////////////////////////////////////////////////////////"
	print tnclone_logo
	print "///////////////////////////////////////////////////////////"
	print 
	print "                                     Author : %s"%(AUTHOR)
	print "                           E-mail : %s\n"%(AUTHOR_EMAIL)
	
def show_configurations(dic):
	
	print "[Options entered for TnClone Analysis]"
	print 
	print "Section 0. Analysis steps selected.\n"
	if dic["sort_on"] :
		print "Sort step       : ON"
	else:
		print "Sort step       : OFF"
		
	if dic["trim_on"] :
		print "Trim step       : ON"
	else:
		print "Trim step       : OFF"
		
	if dic["assembly_on"] :
		print "Assembly step   : ON"
	else:
		print "Assembly step   : OFF"
		
	if dic["analysis_on"] :
		print "Analysis step   : ON"
	else:
		print "Analysis step   : OFF"
	print "\nSection 1. Assembly Options\n"
	print "K-mer size      : %d"%(dic["ksize"])
	if dic["source_seq"]:
		print "Start sequence  : %s"%(dic["source_seq"])
	else:
		print "Start sequence  : NOT SPECIFIED (Maybe obtained from BED file)" 
	if dic["sink_seq"]:
		print "End sequence    : %s"%(dic["sink_seq"])
	else:
		print "End sequence    : NOT SPECIFIED (Maybe obtained from BED file)"
	
	print "Min K-occ       : %d"%(dic["mink_occ"])
	
	if dic["bed_path"] :
		print "Bed file        : %s"%(dic["bed_path"])
		
	else:
		print "Bed file        : NOT SPECIFIED. Start/End is specified"
	
	if dic["down_sample"]:
		print "Down sampling   : YES"
		print "Rate            : %.3f"%(dic["down_sample_rate"])
	else:
		print "Down sampling   : NO"
		print "Rate            : Nan"
		
	if dic["mismatch_correction"]:
		print "Seed correction : MISMATCH"
		print "Num mis allowed : %d"%(dic["max_mis"])
	
	if dic["depth_correction"]:
		print "Seed correction : DEPTH"
		
	if dic["gpu_avail"]:
		print "GPU avail       : YES"
		
	print "\nSection 2. Alignment options\n"
	print "Match score     : %d"%(dic["match_score"])
	print "Mismatch score  : %d"%(dic["mismatch_score"])
	print "Gap open penalty: %d"%(dic["gap_open"])
	print "Gap ext penalty : %d"%(dic["gap_extension"])
	print "Reference path  : %s"%(dic["ref_path"])
	
	
	
	print "\nSection 4. Miscellaneous\n"
	try:
		print "Raw NGS Path    : %s"%(dic["raw_path"])
	except:
		print "Raw NGS Path    : Not specified"
	try:
		print "Sort file       : %s"%(dic["sort_info"])
	except:
		print "Sort file       : Not specified"
	try:
		print "Sample file     : %s"%(dic["sample_info"])
	except:
		print "Sample file     : Not specified"
	try:
		print "Output path     : %s"%(dic["out_path"])
	except:
		print "Output path     : Not specified"
	if dic["draw_coverage"]:
		print "Draw cov. plot  : YES"
	else:
		print "Draw cov. plot  : NO"
		
	if dic["auto"] :
		print "Auto anal. step : YES\n"
	else:
		print "Auto anal. step : NO\n"
def time_format(T):
	hrs = T / 3600
	r1 = T % 3600
	minute = r1 / 60
	sec = r1 % 60
	string = "%d hrs %d min %d sec"%(hrs,minute, sec)
	return string
	
def abspath(path):
	return os.path.abspath(path)

def call_script(script):
	subprocess.call(script ,shell=True)
def launch(script):
	#commands.getoutput(script)
	subprocess.call(script, shell=True)
	
	
class AnalysisThread():
	def __init__(self, opts):
		self.options = opts.option
		self.prefix = "TnClone"
		self.message = []
	def _run_sorter(self):
		outpath = abspath(self.options["out_path"])
		sort_path = outpath + "/sorted"
		
		if not os.path.isdir(sort_path) :
			print "[TnClone::Sorter] Sorting file does not exist. Creating...."
			#logging.info("[LOGGER] Sorting file does not exist. Creating....")
			os.mkdir(sort_path)

		if self.options["sort_on"] == True:
			sort_indir = os.path.abspath(self.options["raw_path"])
			sort_conf = abspath(self.options["sort_info"])
			print "[TnClone::Sorter] Sorting your samples"
			sorter = Sorter(sort_indir,sort_path, sort_conf, self.message)
			#logging.info("[LOGGER] Sorting your samples")
			sorter.sort()
			print "\n".join(self.message)
			
	def _run_trimmer(self):
		outpath = abspath(self.options["out_path"])
		sort_path = outpath + "/sorted"
		trim_path = outpath + "/trimmed"

		if not os.path.isdir(trim_path):
			print "[TnClone::Trimmer] Trimming file path does not exist. Creating...."
			#logging.info("[LOGGER] Trimming file path does not exist. Creating....")
			os.mkdir(trim_path)

		trim_conf = abspath(self.options["sample_info"])
		P = os.path.dirname(os.path.abspath(__file__))
		trimmo_loc = abspath(P + "/Trimmomatic-0.36/trimmomatic-0.36.jar")
		self.message = []
		trimmer = Trimmer(sort_path, trim_path, trim_conf, trimmo_loc,self.message,is_cmd=True)
		print "[TnClone::Trimmer] Trimming your samples"
		trimmer.trim()
		print "\n".join(self.message)
		#logging.info("[LOGGER] Trimming your samples")
		
	def _run_config_gen(self):
		
		logger = []

		outpath = self.options["out_path"]
		assem_conf_name = outpath + "/assembly.config"
		
		if os.path.exists(assem_conf_name) : 
			print "[TnClone::ConfigGenerator] Eventhough config file exist, program will re-create\n"
			os.remove(assem_conf_name)
			#logging.info("[LOGGER] Configuration file already exists. Program will not create.\n")
		else:
			print "[TnClone::ConfigGenerator] Configuration file does not exist (AUTOMODE)\n"
			print "[TnClone::ConfigGenerator] Building assembly configuration file"
		#logging.info("[LOGGER] Configuration file dows not exist (AUTOMODE)\n")
		#logging.info("[LOGGER] Configuration file dows not exist (AUTOMODE)\n")
		if self.options["from_bed"]:
			for line in open(self.options["sample_info"]):
				data = line.rstrip().split()
				gene = data[0]
				nsamples = int(data[1])
				bed_path = self.options["bed_path"]
				ref_path = self.options["ref_path"]
				ksize = int(self.options["ksize"])
				cov = int(self.options["mink_occ"])
				#print gene , nsamples, bed_path, ref_path, ksize, cov
				script = "python " + src_path + "/buildConfigFile.py -k %d -b %s -c %d -n %d -r %s -o %s"%(ksize, bed_path, cov, nsamples, ref_path + "/" + gene + ".gene.fa", assem_conf_name)
				print "[TnClone::ConfigGenerator] " + script
				#logging.info("[LOGGER] " + script)
				launch(script)

		else:
			src = self.options["source_seq"]
			snk = self.options["sink_seq"]
			ksize = int(self.options["ksize"])
			cov = int(self.options["mink_occ"])
			F = open(assem_conf_name , "w")
			for record in open(self.options["sample_info"]):
				data = record.rstrip().split()
				gene = data[0]
				N = int(data[1])
				for x in xrange(N) :
					dat = "%s\t%d\t%s\t%s\t%d\n"%(gene + "_" + str(x+1), ksize, src, snk, cov)
					F.write(dat)

			F.close()
			
	def _run_until_assembly(self):
		outpath = self.options["out_path"]
		
		cov = int(self.options["mink_occ"])

		prefix = self.prefix

		KSIZE = self.options["ksize"]
		ksize = KSIZE

		assem_conf_name = abspath(self.options["out_path"]) + "/assembly.config"

		log = []
		T1 = time.time()
		pid = os.getpid()
		manager = MemoryUsageManager(pid,unit="Gb")
		
		
		with open(assem_conf_name) as conf_file:
			for line in conf_file:
				data = line.rstrip().split("\t")
				try:
					assert len(data) == 5
				except AssertionError:
					raise ValueError("Sth's wring with assem config. E-mail to developer")

				sample = data[0]
				ksize = int(data[1])
				KSIZE = ksize
	
				src = data[2]
				snk = data[3]

				depth_cutoff = data[4]
				gene = sample.split("_")[0]

				
				
				fq_path = outpath + "/trimmed"
				if not os.path.isdir(fq_path):
					print "[TnClone::ERRRR] There is no such a folder(%s) required for assembly. Please see manual throughly."%(fq_path)
					print "[TnClone::ERROR] Program will terminate."
					exit(-1)
				fwd_fqname = fq_path + "/" + gene + "/" + sample + "_R1.trimmed.fastq"
				rev_fqname = fq_path + "/" + gene + "/" + sample + "_R2.trimmed.fastq"
				
				if not os.path.isfile(fwd_fqname):
					print "[TnClone::ERROR] Thene is no such a file(%s) required for assembly. Please see manula throughly"%(fwd_fqname)
					print "[TnClone::ERROR] Program will terminate."
					exit(-1)
					
				if not os.path.isfile(rev_fqname):
					print "[TnClone::ERROR] Thene is no such a file(%s) required for assembly. Please see manula throughly"%(rev_fqname)
					print "[TnClone::ERROR] Program will terminate."
					exit(-1)
				fwd_size = os.stat(fwd_fqname).st_size
				rev_size = os.stat(rev_fqname).st_size

				down_sample = self.options["down_sample"]
				sampling_factor = float(self.options["down_sample_rate"])
				sampling = False
				if down_sample:
					if fwd_size / ( 1024 * 1024 ) >= 30 and rev_size / ( 1024 * 1024 ) >= 30 : ### Greater than 30Mb for uncompressed data
						sampling = True
						print "[TnClone::AssemblyMachinery] Detected file size > 30mb... Downsampling process triggered...\n"
						print "[TnClone::AssemblyMachinery] File size will be reduced from [FWD] %d Mb -> %d Mb || [REV] %d Mb -> %d Mb using factor %f\n"%(fwd_size / ( 1024 * 1024 ) , int(fwd_size / ( 1024 * 1024 ) * sampling_factor) , rev_size / ( 1024 * 1024 ) , int(rev_size / ( 1024 * 1024 ) * sampling_factor) , sampling_factor)
						#logging.info("[LOGGER] Detected file size > 30mb... Downsampling process triggered...\n")
						#logging.info("[LOGGER] File size will be reduced from [FWD] %d Mb -> %d Mb || [REV] %d Mb -> %d Mb using factor %f\n"%(fwd_size / ( 1024 * 1024 ) , int(fwd_size / ( 1024 * 1024 ) * sampling_factor) , rev_size / ( 1024 * 1024 ) , int(rev_size / ( 1024 * 1024 ) * sampling_factor) , sampling_factor))
						sampler = Sampler(fwd_fqname,rev_fqname , sample_ratio=sampling_factor)
						sampler.sample()
					else:
						self.message.append("[TnClone::AssemblyMachinery] Eventhough down-sample is set, file size is smaller than 30 Mb. \nThis will not affect program performance thus there is no down sample process\n")
						#logging.info("[LOGGER] Eventhough down-sample is set, file size is smaller than 30 Mb. \nThis will not affect program performance thus there is no down sample process\n")
						sampling = False
				if sampling:
					fwd_fqname = fq_path + "/" + gene + "/" + sample + "_R1.trimmed.fastq.dwn.fastq"
					rev_fqname = fq_path + "/" + gene + "/" + sample + "_R1.trimmed.fastq.dwn.fastq"

				print "[TnClone::AssemblyMachinery] Building KFS table....\n"
				#logging.info("[LOGGER] Building KFS table....\n")

				kfs_path = outpath + "/kfs"
				if not os.path.isdir(kfs_path):
					print "[TnClone::AssemblyMachinery] KFS does not exist. Create."
					#logging.info("[LOGGER] KFS does not exist. Create.")
					os.mkdir(kfs_path)

				kfs_oname = kfs_path + "/" + prefix + "_" + sample + "_k" + str(ksize) + ".kfs"

				kfs = KFS(kfs_oname, ksize=ksize)	

				srp = SRP(fwd_fqname, rev_fqname)

				srp.readFwd()
				srp.readRev()

				for fwd_read in srp.fwdread:
					kfs.computeFreq(fwd_read)

				for rev_read in srp.revread:
					kfs.computeFreq(rev_read)
		
				kfs.writeRecord()

				del srp
				del kfs


				print "[TnClone::AssemblyMachinery] Build Graph using Illumina sequencing reads"
				#logging.info("[LOGGER] Build Graph using Illumina sequencing reads")
				tmp_path = outpath + "/tmp"

				if not os.path.isdir(tmp_path):
					print "[TnClone::AssemblyMachinery] temporal directory not exist. Create"
					#logging.info("[LOGGER] temporal directory not exist. Create")
					os.mkdir(tmp_path)

				assem_path = outpath + "/assem"
				if not os.path.isdir(assem_path):
					print "[TnClone::AssemblyMachinery] No assembly result directory. Create..."
					#logging.info("[LOGGER] No assembly result directory. Create...")
					os.mkdir(assem_path)
				cfastq_path = outpath + "/cfastq"

				if not os.path.isdir(cfastq_path):
					print "[TnClone::AssemblyMachinery] No CFASTQ directory. Create"
					#logging.info("[LOGGER] No CFASTQ directory. Create")
					os.mkdir(cfastq_path)

				node_tmp = tmp_path + "/" + prefix + "_" + sample + "_k" + str(ksize) + ".nodes"

				assem_outname = assem_path + "/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig"
				mink_occ = int(self.options["mink_occ"])
				assembler = Assembler(kfs_oname, src, snk,  assem_outname, ksize, fwd_fqname, rev_fqname, node_tmp, mink_occ)

				print "[TnClone::AssemblyMachinery] Building...."
				print "[TnClone::AssemblyMachinery] Building graph...\n"
				#sys.stderr.write("[MACHINERY] Building graph...\n")
				#logging.info("[LOGGER] Building....")
				#logging.info("[MACHINERY] Building graph...\n")
				
				assembler.buildGraph()
				
				print "[TnClone::AssemblyMachinery] Dump raw connection information...."
				#logging.info("[LOGGER] Dump raw connection information....")
				#log.append("[LOGGER] Dump connection information....")
				assembler.dump2Tmp(node_tmp)
				print "[TnClone::AssemblyMachinery] graph modification...."
				#logging.info("[LOGGER] graph modification....")
				assembler.modifyGraph(node_tmp)
				print "[TnClone::AssemblyMachinery] Dump connection information...."
				#logging.info("[LOGGER] Dump connection information....")
				assembler.dump2Tmp(node_tmp)
				
				
				print "[TnClone::AssemblyMachinery] Launching Assembly"
				print "[TnClone::AssemblyMachinery] Assembly starts...\n"
				
				#logging.info("[LOGGER] Launching Assembly")
				#logging.info("[MACHINERY] Assembly starts...\n")
				#log.append("[LOGGER] Launching Assembly")
				#sys.stderr.write("[MACHINERY] Assembly starts...\n")
				#sys.stderr.flush()
				assem_out = open(assem_outname, "w")
				ncontig = 0
				assembly_seed_error=False
				### before assembly, check for proper seed
				try:
					assembler.graph.bdg[src]

				except KeyError:
					assembly_seed_error = True
					G = assembler.graph.bdg
					seed = src
					path = outpath
					err_log = self.message
					max_mis = self.options["max_mis"]
					use_depth_correction = self.options["depth_correction"]
					mismatch_correction = self.options["mismatch_correction"]
					gpu_avail = self.options["gpu_avail"]

					ESC = ErrorSeedCorrector(G , seed, path, err_log, max_mis=max_mis,use_depth_profile=use_depth_correction, nvidia_enable = gpu_avail)
					ESC.select()

				if not assembly_seed_error:	
					try:
						
						for path in assembler.assemble().walk():
							#log.append(path)
							if "NOSEED" in path:
								assem_out.write("NO SEED CONTIG\n")
								continue
							if "BROKEN" in path:
								niter = path[-1]
								total_depth , seq ,node_used, d_list , r_dict = path[0], path[1][0], path[2] , path[3], path[4]
								for entry in path[1][1:]:
									kmers = entry
									seq += kmers[-1]
								if seq != '':
									assem_out.write(">BROKEN PATH|%d-th iteration\n"%(niter+1))
									assem_out.write(seq + "\n")

								else:				
									assem_out.write(">BROKEN PATH|%d-th iteration\n"%(niter+1))
								continue

							total_depth , seq ,node_used, d_list , r_dict = path[0], path[1][0], path[2] , path[3], path[4]
							
							r_list = r_dict.values()
							#print r_list
							
							for entry in path[1][1:] :
								kmers = entry
								seq += kmers[-1]

							if seq != '' :
								### 2017-03-16 modified 
								
								Statistics = ContigStats2(seq , int(self.options["ksize"]) , assembler.graph.bdg)
								path_score = Statistics.path_score()
								avg_depth , cv = Statistics.cv()
								assem_out.write(">contig:%d|avg_dp:%d|path_score:%f|cv:%f|sample_name:%s\n"%(ncontig, avg_depth,path_score, cv, sample))
								assem_out.write(seq + "\n")
								
								ncontig += 1
								
					except IOError:
						assembler.clearAssembler()
						continue
					except KeyError:
						assembler.clearAssembler()
						#assem_out.write("NO SEED CONTIG\n")
					assem_out.close()
					assembler.clearAssembler()
					print "[TnClone::AssemblyMachinery] Assembly done...\n"
					#logging.info("[MACHINERY] Assembly done...\n")

				else:
					check = False
					if os.stat(outpath + "/seed.fa") == 0 :
						assem_out.write("NO SEED CONTIG\n")
						continue
					for line in open(outpath + "/seed.fa") :
						src = line.rstrip()
						assembler.set_src(src)### change source sequence for assembly
						try:
							for path in assembler.assemble().walk():
								#log.append(path)
								if "NOSEED" in path:
									assem_out.write("NO SEED CONTIG\n")
									continue
								if "BROKEN" in path:
									niter = path[-1]
									total_depth , seq ,node_used, d_list , r_dict = path[0], path[1][0], path[2] , path[3], path[4]
									for entry in path[1][1:]:
										kmers = entry
										seq += kmers[-1]
									if seq != '':
										assem_out.write(">BROKEN PATH|%d-th iteration\n"%(niter+1))
										assem_out.write(seq + "\n")
	
									else:				
										assem_out.write(">BROKEN PATH|%d-th iteration\n"%(niter+1))
									continue

								total_depth , seq ,node_used, d_list , r_list = path[0], path[1][0], path[2] , path[3], path[4]
								#print r_list
								
								for entry in path[1][1:] :
									kmers = entry
									seq += kmers[-1]

								if seq != '' :
									Statistics = ContigStats2(seq , int(self.options["ksize"]) , assembler.graph.bdg)
									path_score = Statistics.path_score()
									avg_depth , cv = Statistics.cv()
									assem_out.write(">contig:%d|avg_dp:%d|path_score:%f|cv:%f|sample_name:%s\n"%(ncontig, avg_depth,path_score, cv, sample))
									assem_out.write(seq + "\n")
									ncontig += 1
									check = True
						except IOError:
							assembler.clearAssembler()
							continue
						except KeyError:
							assembler.clearAssembler()
							#assem_out.write("NO SEED CONTIG\n")
							continue
						if check == True :
							break ## excape for loop
					assem_out.close()
					assembler.clearAssembler()
					print "[TnClone::AssemblyMachinery] Assembly done...\n"
					#logging.info("[MACHINERY] Assembly done...\n")
				#sys.stderr.write("[MACHINERY] Assembly done...\n")
				#sys.stderr.flush()

				### REQUIRED : DISCARD CONTIGS
			### sixth process : Downstream analysis
				ss_find_file = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.ss.rmv"
				ssf = SSF(assem_outname , ss_find_file, src, snk,self.message)
				ssf.run()

				dup_rmv_file = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.final.fa"
				drt = DRT(ss_find_file , dup_rmv_file,self.message)
				drt.run()
		T2 = time.time()
		manager.profile()
		print "[TnClone::AssemblyMachinery] MEM used : " + manager.mem_usage + "\n"
		print "[TnClone::AssemblyMachinery] Performed in %f sec\n"%(T2-T1)
		
	def _run_alignment(self):
		m = int(self.options["match_score"])
		x = int(self.options["mismatch_score"])
		g = int(self.options["gap_open"])
		e = int(self.options["gap_extension"])

		outpath = abspath(self.options["out_path"])
		refpath = abspath(self.options["ref_path"])

		log = []

		sam_path = outpath + "/sam"
		prefix = self.prefix
		if not os.path.isdir(sam_path) :
			print "[TnClone::Aligner] SAM folder do not exist. Create."
			#logging.info("[LOGGER] SAM folder do not exist. Create.")
			os.mkdir(sam_path)
		outpath = abspath(self.options["out_path"])
		vcf_path = outpath + "/vcf"
		#log = []
		prefix = self.prefix
		if not os.path.isdir(vcf_path):
			print "[TnClone::Caller] No VCF dir. Create."
			#logging.info("[LOGGER] No VCF dir. Create.")
			os.mkdir(vcf_path)

		assem_fname = abspath(outpath + "/assembly.config")
		no_seed_samples = []
		broken_samples = []
		with open(assem_fname) as conf_file:
			for line in conf_file:
				data = line.rstrip().split("\t")

				try:
					sample = data[0]
					ksize = int(data[1])

					refname = refpath + "/" + sample.split("_")[0] + ".gene.fa"

					query = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.final.fa"

					tmp_f = open(query,"r")

					entry = tmp_f.readlines()

					if len(entry) == 0 : continue

					if entry[0].rstrip() == "NO SEED CONTIG" :
						no_seed_samples.append(sample)
						tmp_f.close()

						continue

					if entry[0].rstrip().startswith(">BROKEN PATH"):
						broken_samples.append(sample)
						tmp_f.close()
						continue

					samname = sam_path + "/" + sample + "_k" + str(ksize) + ".dna.sam"

					script = src_path + "/SSW-CPP -m %d -x %d -g %d -e %d %s %s > %s"%(m,x,g,e,refname,query,samname)

					print "[TnClone::Aligner] " + script+"\n"
					#logging.info(script+"\n")

					launch(script)
					vcfname = vcf_path + "/" + sample + "_k" + str(ksize) + ".dna.vcf"
					script2 = "python " + src_path + "/variant_caller.py --sam %s --ref %s --vcf %s"%(samname, refname, vcfname)
					print "[TnClone::Caller] " + script2 + "\n"
					#logging.info(script2 + "\n")
					launch(script2)
					aa_contigname = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.final.fa.aa.fa"
					aa_refname = refpath + "/" + prefix + "_" + sample.split("_")[0] + ".gene.fa.aa.fa"

					aa_contigf = open(aa_contigname , "w")
					aa_reff = open(aa_refname , "w")


					if seqio_avail :
						for record in SeqIO.parse(refname , "fasta") :
							id = record.id
							seq = str(record.seq)
							### Considering frame
							#for frame in [0,1,2]:
							#	acid = str(Seq(seq[frame:]).translate())
							acid = str(Seq(seq).translate())
							aa_reff.write(">" + "\n" + acid + "\n")
						aa_reff.close()

						for record in SeqIO.parse(query,"fasta") :
							id = record.id
							seq = str(record.seq)
							for frame in [0,1,2]:
								acid = str(Seq(seq[frame:]).translate())
								aa_contigf.write(">" + id + "-frame:" + str(frame) + "\n" + acid + "\n")
						aa_contigf.close()

					else:
						ref_parser = FastaParser(refname)
						contig_parser = FastaParser(query)
						ref_parser.open()
						ref_parser.read()
						rid = ref_parser.id
						aa = translate(ref_parser.seq)
						aa_reff.write(rid + "\n" + aa + "\n")
						aa_reff.close()
						ref_parser.close()
				
						contig_parser.open()
						while 1:
							try:
								contig_parser.read()
								cid = contig_parser.id
								for frame in [0,1,2]:
									caa = translate(contig_parser.seq[frame:])
									aa_contigf.write(cid + "-frame:" + str(frame) + "\n" + caa + "\n")
							except StopIteration:
								break
						aa_contigf.close()		
						contig_parser.close()
					acid_samname = sam_path + "/" + sample + "_k" + str(ksize) + ".protein.sam"
					script3 = src_path + "/SSW-C -m %d -x %d -o %d -e %d -p -a %s/blosum62.txt -c -s -h %s %s > %s"%(m,x,g,e,src_path,aa_refname , aa_contigname, acid_samname)
					print "[TnClone::Aligner] " + script3 + "\n"
					#logging.info(script3 + "\n")
					launch(script3)
					aa_vcfname = vcf_path + "/" + sample + "_k" + str(ksize) + ".protein.vcf"
					script4 = "python " + src_path + "/variant_caller.py --sam %s --ref %s --vcf %s"%(acid_samname, aa_refname, aa_vcfname)
					print "[TnClone::Caller] " + script4 + "\n"
					#logging.info(script4 + "\n")
					launch(script4)
				except IOError:
					continue
		sam_path = outpath + "/sam"

		pass_dict = {}
		with open(self.options["sample_info"]) as info_file:
			for line in info_file :
				data = line.rstrip().split("\t")
				pass_dict[data[0]] = int(data[1])

		assem_path = outpath + "/assem"
		
		report_path = outpath + "/report"
		if not os.path.isdir(report_path):
			print "[TnClone::Reporter] No Report dir.Create."
			#logging.info("[LOGGER] No Report dir.Create.")
			os.mkdir(report_path)

		report_name = outpath + "/report/analysis.report"

		now = time.localtime()
		analysis_time = "%04d-%02d-%02d %02d:%02d:%02d" % (now.tm_year, now.tm_mon, now.tm_mday, now.tm_hour, now.tm_min, now.tm_sec)
		report = open(report_name , "aw")
		authorship = """############################################
ANALYSIS INFORMATION
CREATED BY TnClone
PROGRAM AUTHOR : SUNGHOON HEO
PLEASE CONTACT AUTHOR VIA E-MAIL BELOW IF YOU HAVE ISSUE
tahuh1124(at)gmail.com
###########################################\n"""
		report.write(authorship)
		report.write("ANALYSIS FINISHED AT : %s\n"%(analysis_time))
		report.close()

		CNR = ContigNumberReporter(assem_path, report_name, pass_dict, prefix)
		KSIZE = int(self.options["ksize"])
		CNR.run(KSIZE)
		"""
		AAR = AminoAcidReporter(sam_path, report_name, pass_dict, prefix,no_seed_samples,broken_samples)

		AAR.init()
		AAR.analyse(KSIZE)

		DNAR = DNAReporter(sam_path, report_name, pass_dict , prefix,no_seed_samples,broken_samples)

		DNAR.init()

		DNAR.analyse(KSIZE)
		"""
		dna_reporter = ErrorReporter(outpath, prefix, KSIZE , type='dna')
		protein_reporter = ErrorReporter(outpath, prefix, KSIZE , type='protein')
		with open(self.options["sample_info"]) as info_file:
			for line in info_file:
				data = line.rstrip().split("\t")
				gene = data[0]; nsamples = int(data[1])
				for x in xrange(1,nsamples + 1):
					dna_reporter.analyse(gene , str(x))
					protein_reporter.analyse(gene, str(x))

			dna_reporter.write_result()
			protein_reporter.write_result()
		
		
		
		
		#logging.info("#####################################")
		#logging.info("[TnClone] Finished whole analysis. Please enjoy result!")
		#logging.info("#####################################")
		
	
		
		
class DiagnosisThread:
	def __init__(self, sample_info_file, trim_path, ref_path):
		self.sample_info = sample_info_file
		self.trim_path = trim_path
		self.ref_path = ref_path
		self.log_prefix = "[TnClone::DiagnosisPlot] "
		self.root_path = "/".join(trim_path.split("/")[:-1])
		
	def run_analysis(self):
		diagnosis_path = self.root_path + "/diagnosis"
		sam_diagnosis = diagnosis_path + "/sam"
		bam_diagnosis = diagnosis_path + "/bam"
		pileup_diagnosis = diagnosis_path + "/pileup"
		if not os.path.isdir(diagnosis_path):
			print self.log_prefix + "Diagnosis folder do not exist. Create."
			os.mkdir(diagnosis_path)

		if not os.path.isdir(sam_diagnosis):
			print self.log_prefix + "Diagnosis->SAM folder do not exist. Create."
			os.mkdir(sam_diagnosis)

		if not os.path.isdir(bam_diagnosis):
			print self.log_prefix + "Diagnosis->BAM folder do not exist. Create."
			os.mkdir(bam_diagnosis)

		if not os.path.isdir(pileup_diagnosis):
			print self.log_prefix + "Diagnosis->Pileup folder do not exist. Create."
			os.mkdir(pileup_diagnosis)
		cnt = 0
		with open(self.sample_info) as F:
			for line in F:
				data = line.rstrip().split()
				n = int(data[1])
				cnt += 1

		print self.log_prefix + "Number of samples are {}. Diagnosis will take a while".format(cnt)
		t1 = time.time()
		with open(self.sample_info) as F :
			for line in F :
				data = line.rstrip().split("\t")
				gene = data[0]
				numbers = int(data[1])
				REF = self.ref_path + "/" + gene + ".gene.fa"
				bwa_index = src_path + "/bwa index %s"%(REF)
				print self.log_prefix + bwa_index
				call_script(bwa_index)

				for i in xrange(1,numbers + 1):
					### check for file existance
					### Check if down streamed
					if os.path.isfile(self.trim_path + "/" + gene + "/" +gene + "_" + str(i) + "_R1.trimmed.fastq.dwn.fastq"):
						fwd_fastq = self.trim_path + "/" + gene + "/" +gene + "_" + str(i) + "_R1.trimmed.fastq.dwn.fastq"
						rev_fastq = self.trim_path + "/" + gene + "/" +gene + "_" + str(i) + "_R2.trimmed.fastq.dwn.fastq"
					else:
						fwd_fastq = self.trim_path + "/" + gene + "/" + gene + "_" + str(i) + "_R1.trimmed.fastq"			
						rev_fastq = self.trim_path + "/" + gene + "/" + gene + "_" + str(i) + "_R2.trimmed.fastq"
					bwa_sam = sam_diagnosis + "/" + gene + "_" + str(i) + ".bwa.sam"
					bwa_mem = src_path + "/bwa mem -M -R '@RG\tID:%s\tSM:%s\tLB:%s\tPL:ILLUMINA\tPU:%s' %s %s %s > %s"%(gene,gene,gene,gene,REF,fwd_fastq,rev_fastq,bwa_sam)
					print self.log_prefix + bwa_mem
					call_script(bwa_mem)
					bam1 = bam_diagnosis + "/" + gene + "_" + str(i) + ".bam"
					samtools_view = src_path + "/samtools view -b -T %s -S %s -o %s"%(REF,bwa_sam,bam1)
					print self.log_prefix + samtools_view
					call_script(samtools_view)
					bam2 = bam1 + ".sorted"
					samtools_sort = src_path + "/samtools sort %s %s"%(bam1,bam2)
					print self.log_prefix + samtools_sort
					call_script(samtools_sort)
					samtools_index = src_path + "/samtools index %s %s"%(bam2+ ".bam", bam2 + ".bam.bam")
					print self.log_prefix + samtools_index
					call_script(samtools_index)
					
					pileup_file = pileup_diagnosis + "/" + gene + "_" + str(i) + ".pileup"
					samtools_pileup = src_path + "/samtools mpileup -x -q 1 -Q 1 -d 1000000 -s -f %s %s > %s"%(REF, bam2+".bam",pileup_file)
					print self.log_prefix + samtools_pileup
					call_script(samtools_pileup)
					t2 = time.time()
		
		t2 = time.time()
		print self.log_prefix + "Done diagnosis. Elapsed {} seconds".format(t2-t1)
"""
Compared to TnClone's GUI version, CMD version will draw all plots available
"""
class DiagnosisPlot:
	def __init__(self, sample_info_file, trim_path, ref_path):
		self.sample_info = sample_info_file
		self.trim_path = trim_path
		self.ref_path = ref_path
		self.log_prefix = "[TnClone::DiagnosisPlot] "
		self.root_path = "/".join(trim_path.split("/")[:-1])
		self.dataset = None
		self.first_draw = False
		self.filenames = None
		self.did_diagnosis = False
		self.log_screen = []
	def diagnosis(self):
		self.analyser = DiagnosisThread(self.sample_info ,self.trim_path, self.ref_path)
		
		if os.path.isdir("/".join(self.trim_path.split("/")[:-1]) + "/diagnosis/pileup"):
			self.did_diagnosis = True
		if not self.did_diagnosis :
			self.analyser.run_analysis()
	def _parse_all_pileups(self) :
		### return all parsed data
		pileup_path = "/".join(self.trim_path.split("/")[:-1]) + "/diagnosis/pileup"
		files = os.listdir(pileup_path)
		dataset = {}
		for filename in files :
			print self.log_prefix + " Pileup-file : %s"%(filename)
			dataset[filename] = {"posits" : [] , "depths": []}

			with open(pileup_path + "/" + filename) as PILEUP:
				for line in PILEUP:
					### storing data
					data = line.rstrip().split("\t")
					position = int(data[1])
					depth = int(data[3])
					dataset[filename]["posits"].append(position)
					dataset[filename]["depths"].append(depth)
		self.dataset = dataset


	def draw(self):
		### Diagnosis first
		
		self.diagnosis()
			
		self._parse_all_pileups()
		
		dataset = self.dataset
		filenames = dataset.keys()
		self.filenames = filenames
		
		for fname in self.filenames :
		
			data = dataset[fname]
			
			xmin = min(data['posits'])
			xmax = max(data['posits'])
			
			fig = plt.figure(figsize=(10,8.5))
			
			subplt = fig.add_subplot(111)
			subplt.set_ylabel('Coverage')
			subplt.set_xlabel('Coordinates')
			axes = plt.gca()
			axes.set_xlim([xmin,xmax])
			posits = data['posits']
			depths = data['depths']
			
			subplt.plot(posits, depths , lw = 2)
			
			subplt.axhline(y=np.mean(depths) , ls = 'dashed' , color = 'red' , lw = 3)
			subplt.fill_between(posits, depths, facecolor='blue', alpha=0.5)


			for label in subplt.get_yticklabels():
				label.set_visible=(False)
			fig.suptitle('Diagnostic plot for evaluating coverage of given reference\n(Dashed line indicates average coverage)')
			
			save_fig_name = self.root_path + "/diagnosis/" + fname + '-coverage.png'
			
			plt.savefig(save_fig_name)


		
class Option(object):
	def __init__(self, option_dict=None):
		self.option = option_dict
		
	
if __name__ == "__main__" :
	show_header_program()
	
	parser = optparse.OptionParser(description='TnClone CUI version')
	parser.add_option("--ngs_path" , \
						default = None,
						type=str, \
						help = "Raw NGS result path")
						
	parser.add_option("--sort_file", \
						type=str,\
						default = None,
						help="Sorting information containing file.")
	parser.add_option("--sample_file", \
						type=str,\
						default = None,
						help = "Sample information file.")
						
	parser.add_option("--start",\
						type=str,\
						default = None,
						help="Start sequence of assembly. Length must be same as k-mer size option")
						
	parser.add_option("--end",\
						type=str,\
						default = None,
						help="End sequence of assembly. Length must be same as k-mer size option")
						
	parser.add_option("--ref_path",\
						type=str,\
						default = None,
						help="Reference directory that contains reference files for analysis")
						
	parser.add_option("--output_path", \
						type=str,\
						default = None,
						help = "Output path")
						
	parser.add_option("--kmer_size" ,\
						type=int,\
						help='K-mer size. Defaul = 63',\
						default=63)
						
	parser.add_option("--min_kocc",\
						type=int,\
						help="Minimum k-mer occurence. \
						This value is used to reduce amboguity. Default = 3",\
						default=3)
						
	parser.add_option("--sort_on", \
						type=int,\
						help="Check whether to turn on sorting. o for no 1 for yes. default = 1",
						default = 1)
	parser.add_option("--trim_on", \
						type=int,\
						help="Check whether to turn on trimming. o for no 1 for yes. default = 1",
						default = 1)

	parser.add_option("--assembly_on", \
						type=int,\
						help="Check whether to turn on assembly. o for no 1 for yes. default = 1",
						default = 1)

	parser.add_option("--analysis_on", \
						type=int,\
						help="Check whether to turn on analysis. o for no 1 for yes. default = 1",
						default = 1)

	parser.add_option("--auto_analysis", \
						type=int,\
						help="Generates analysis automatically. Will start from analysis you specified\
						0 for no 1 for yes. default = 1",\
						default = 1)
	parser.add_option("--bed_file" ,\
						type=str,\
						help="Bed file",
						default=None)
						
	parser.add_option("--down_sampling", \
						type=int,\
						help="Downsampling or not 1 for yes and 0 for not. Default = 0",\
						default=0)
						
	parser.add_option("--down_sample_rate", \
						type=float,\
						help="Down sampling ratio. Default = 0.3",\
						default=0.3)
						
	parser.add_option("--seed_mismatch_correction_method",\
						type=str,\
						help="Type of methode to correct seed error. One of depth or mismatch.\
						default=mismatch",\
						default="mismatch")
						
	parser.add_option("--num_mismatch",\
						type=int,\
						help="Number of mismatch allowed. \
						if seed mismatch correction method is set depth , this will not invoked.\
						Default = 3",\
						default=3)
	parser.add_option("--gpu_avail", \
						type=int,\
						help="Let the program know if gpu (NVidia card) avail. \
						0 for no 1 for yes. Default = 0",\
						default=0)
						
	parser.add_option("--aln_match", \
						type=int,\
						help="Match score. INT value. Default = 1",\
						default=1)
						
	parser.add_option("--aln_mismatch",\
						type=int,\
						help="Mismatch score. INT value. Will work as negative value. Default = 3",\
						default=3)
						
	parser.add_option("--aln_open_penalty" ,\
						type=int,\
						help="Gap open penalty. INT value. Default = 5",\
						default=5)
						
	parser.add_option("--aln_gap_extension", \
						type=int,\
						help="Gap extension penalty. INT value. Default = 3",\
						default=3)
	
	parser.add_option("--draw_coverage", \
						type=int,\
						help="Draw coverage plot before analyis. 0 for no 1 for yes. Default = 0",\
						default=0)
						
	opt , args = parser.parse_args()
	
	ngs_path = opt.ngs_path; sort_file = opt.sort_file; sample_file = opt.sample_file
	start = opt.start; end = opt.end; ref_path = opt.ref_path
	output_path = opt.output_path; kmer_size = opt.kmer_size; min_kocc = opt.min_kocc
	sort_on = opt.sort_on; trim_on = opt.trim_on; assem_on = opt.assembly_on; anal_on = opt.analysis_on
	bed_file = opt.bed_file; down_sampling = opt.down_sampling
	down_sample_rate = opt.down_sample_rate; seed_coor_method = opt.seed_mismatch_correction_method
	num_mismatch = opt.num_mismatch; gpu_avail = opt.gpu_avail
	aln_match = opt.aln_match; aln_mismatch = opt.aln_mismatch
	aln_open_penalty = opt.aln_open_penalty; aln_ext_penalty = opt.aln_gap_extension
	
	draw_coverage = opt.draw_coverage
	
	auto_anal_gen = opt.auto_analysis
	
	if len(sys.argv) < 2 :
		#print "[TnClone::Main] There is not options tou have passed.\n\
		#python tnclone-cmd.py -h to see options."
		show_help()
		print "\nNo options are specified. Program will terminate."
		exit(-1)
		
	### Configuring
	
	### Some Flag values
	raw_exist = False
	sort_file_exist = False
	sample_file_exist = False
	source_exist = False
	sink_exist = False
	bed_exist = False
	ref_exist = False
	output_exist = False
	kmer_len_same = False
	analysis_mark = False
	down_sample_rate_set = False
	seed_coor_mismatch_set = False
	seed_coor_mismatch_value_set = False
	### Read dict
	opt_dict ={}
	if ngs_path :
		opt_dict["raw_path"] = ngs_path
		raw_exist = True
	if sort_file :
		opt_dict["sort_info"] = sort_file
		sort_file_exist = True
	if sample_file:
		opt_dict["sample_info"] = sample_file
		sample_file_exist = True
	if start :
		opt_dict["source_seq"] = start
		source_exist = True
	else:
		opt_dict["source_seq"] = None
	if end:
		opt_dict["sink_seq"] = end
		sink_exist = True
	else:
		opt_dict["sink_seq"] = None
	if ref_path :
		opt_dict["ref_path"] = ref_path
		ref_exist = True
	else:
		opt_dict["ref_path"] = None
	if output_path :
		opt_dict["out_path"] = output_path
		output_exist = True
	else:
		opt_dict["out_path"] = None
	if kmer_size:
		opt_dict["ksize"] = kmer_size
	if min_kocc:
		opt_dict["mink_occ"] = min_kocc
	
	
	if sort_on:
		
		analysis_mark = True
		opt_dict["sort_on"] = True
		if auto_anal_gen:
			opt_dict["trim_on"] = True
			opt_dict["assembly_on"] = True
			opt_dict["analysis_on"] = True
		else:
			opt_dict["trim_on"] = False
			opt_dict["assembly_on"] = False
			opt_dict["analysis_on"] = False
	else:
		opt_dict["sort_on"] = False
	if trim_on:
		
		analysis_mark = True
		opt_dict["trim_on"] = True
		if auto_anal_gen:
			opt_dict["assembly_on"] = True
			opt_dict["analysis_on"] = True
		else:
			opt_dict["assembly_on"] = False
			opt_dict["analysis_on"] = False
	else:
		opt_dict["trim_on"] = False
	if assem_on:
		
		analysis_mark = True
		opt_dict["assembly_on"] = True
		if auto_anal_gen:
			opt_dict["analysis_on"] = True
		else:
			opt_dict["analysis_on"] = False
	else:
		opt_dict["assembly_on"] = False
	if anal_on:
		if not auto_anal_gen:
			opt_dict["analysis_on"] = False
		else:
			analysis_mark = True
			opt_dict["analysis_on"] = True
	else:
		opt_dict["analysis_on"] = False
			
	if bed_file != None :
		opt_dict["bed_path"] = bed_file
		opt_dict["from_bed"] = True
		bed_exist = True
	else:
		opt_dict["bed_path"] = None
		opt_dict["from_bed"] = False
	if down_sampling:
		opt_dict["down_sample"] = True
		if down_sample_rate:
			opt_dict["down_sample_rate"] = down_sample_rate
			down_sample_rate_set = True
		else:
			down_sample_rate_set = False
	else:
		opt_dict["down_sample"] = False	
		opt_dict["down_sample_rate"] = down_sample_rate
	
	if seed_coor_method:
		if seed_coor_method == "mismatch":
			opt_dict["mismatch_correction"] = True
			opt_dict["depth_correction"] = False
			seed_coor_mismatch_set = True
		if seed_coor_method == "depth":
			opt_dict["depth_correction"] = True
			opt_dict["mismatch_correction"] = False
			
	if num_mismatch != None:
		opt_dict["max_mis"] = num_mismatch
		seed_coor_mismatch_value_set = True
		
	if gpu_avail:
		opt_dict["gpu_avail"] = True
	else:
		opt_dict["gpu_avail"] = None
	if aln_match != None:
		opt_dict["match_score"] = aln_match
	if aln_mismatch != None:
		opt_dict["mismatch_score"] = aln_mismatch
	if aln_open_penalty != None:
		opt_dict["gap_open"] = aln_open_penalty
	if aln_ext_penalty != None:
		opt_dict["gap_extension"] = aln_ext_penalty
		
	opt_dict["draw_coverage"] = draw_coverage
	opt_dict["auto"] = auto_anal_gen
	### SHOW configuration
	
	show_configurations(opt_dict)
	
	
	### Check if there is configuration errors
	error_prefix = "[TnClone-CMD::Error] "
	### If source and sink is not set and bed file not set
	if (not start) and (not end):
		if not bed_file:
			print error_prefix + "Start/End sequence must be specified or bed file must be specified"
			print error_prefix + "Program will terminate."
			exit(-1)
	
	### Check if source ans sink length is equal to k-mer size
	if  start:
		if len(start) != kmer_size:
			print error_prefix + "Start sequence length is differ from k-mer length setting"
			print error_prefix + "Start sequence length : k-mer setting = %d : %d"%(len(start), kmer_size)
			print error_prefix + "Program will terminate."
			exit(-1)
	if end:
		if len(end) != kmer_size:
			print error_prefix + "End sequence length is differ from k-mer length setting"
			print error_prefix + "End sequence length : k-mer setting = %d : %d"%(len(end), kmer_size)
			print error_prefix + "Program will terminate."
			exit(-1)
		
	### If analysis option sort is set but no raw ngs path set
	if sort_file_exist:
		if not raw_exist:
			print error_prefix + "To sort your data, raw NGS path must be specified."
			print error_prefix + "Program will terminate."
			exit(-1)
			
	### If other analysis option is set and sample file did not specified
	if opt_dict["trim_on"] :
		if not sample_file_exist:
			print error_prefix + "To trim your data, sample file must be specified"
			print error_prefix + "Also check if your files are already sorted. See Manual for details"
			print error_prefix + "Program will terminate."
			exit(-1)
			
	if opt_dict["assembly_on"] :
		if not sample_file_exist:
			print error_prefix + "To launch assembly, sample file must be specified"
			print error_prefix + "Program will terminate."
			exit(-1)
		
	if opt_dict["analysis_on"]:
		if not sample_file_exist:
			print error_prefix + "To launch down stream analysis, sample file must be specified"
			print error_prefix + "Program will terminate."
			exit(-1)
			
		if not ref_exist:
			print error_prefix + "To launch down stream analysis, reference file must be specified"
			print error_prefix + "Program will terminate."
			exit(-1)
			
	### If output path not specified
	if not output_exist:
		print error_prefix + "Need output path for analysis."
		print error_prefix + "Program will terminate."
		exit(-1)
	
	### When analysis step did not specified
	if not analysis_mark:
		print error_prefix + "Please specify analysis step."
		print error_prefix + "Program will terminate."
		exit(-1)
	
	### Seed method check
	if seed_coor_mismatch_set:
		if not seed_coor_mismatch_value_set:
			print error_prefix + "Seed correction method is set as mismatch \
			but the numbers are not set. Please set numbers"
			print error_prefix + "Program will terminate."
			exit(-1)
			
	### Downsample arguments check
	if opt_dict["down_sample"]:
		if not down_sample_rate_set:
			print error_prefix + "Downsample option is set but the rate is not set. Please specify"
			print error_prefix + "Program will terminate."
			exit(-1)
			
	### GPU check
	if opt_dict["gpu_avail"]:
		nvcc_check = "nvcc 2>./nvcc_check.txt"
		launch(nvcc_check)
		with open("./nvcc_check.txt") as f:
			l = f.readline()
			if l.startswith("The program 'nvcc' is currently not installed"):
				print error_prefix + "Nvida card not detected or CUDA is not installed. Please check"
				print error_prefix + "Program will terminate."
				exit(-1)
				
	if draw_coverage:
		if not sample_file_exist:
			print error_prefix + "Drawing option is set but there is no sample information."
			print error_prefix + "Program will terminate."
			exit(-1)
			
		if not output_exist:
			print error_prefix + "Drawing option is set but there is no output path information"
			print error_prefix + "Program will terminate."
			exit(-1)
			
		if not ref_exist:
			print error_prefix + "Drawing option is set but there is no reference information"
			print error_prefix + "Program will terminate."
			exit(-1)
		if not os.path.isdir(output_path + "/trimmed"):
			print error_prefix + "Before analysis, FASTQ files must be trimmed.\n\
[Solutions]\
\n1. Please first lauch this program using option --trim_on=1 (if you have already sorted, else requires --sort_on=1 option also.). Then re-lauch this program with\
drawing option.\n\
2. Locate your fastq files under output_path/trimmed/[target_gene_name]/ with names \
\nFor first paird end read \
[targe_gene_name]_[sample_number]_R1.trimmed.fastq \
\nFor second paired end read\
[targe_gene_name]_[sample_number]_R2.trimmed.fastq \
\nSee Manual for details."
			print error_prefix + "Program will terminate."
			exit(-1)
		print "[TnClone::DiagnosisPlot] Drawing your coverage images at %s"%(output_path + "/diagnosis/")	
		drawer = DiagnosisPlot(sample_file, output_path + "/trimmed", ref_path)
		drawer.draw()
		
	#### Now it is time to trigger TnClone!!!
	opt = Option(opt_dict)
	pipeline_thread = AnalysisThread(opt)
	
	t1 = time.time()
	if pipeline_thread.options["sort_on"]:
		pipeline_thread._run_sorter()
		
	if pipeline_thread.options["trim_on"]:
		pipeline_thread._run_trimmer()
		
	if pipeline_thread.options["assembly_on"]:
		pipeline_thread._run_config_gen()
		pipeline_thread._run_until_assembly()
		
	if pipeline_thread.options["analysis_on"]:
		pipeline_thread._run_alignment()
		
	t2 = time.time()
	
	elapsed = t2 - t1
	formatted_time_str = time_format(elapsed)
	
	print "[TnClone::Timer] Elapsed %s"%(formatted_time_str)
	
	print "#########################################################"
	print "#                                                       #"
	print "#                                                       #"
	print "#[TnClone] Finished whole analysis. Please enjoy result!# "
	print "#                                                       #"
	print "#                                                       #"
	print "#########################################################"
		
	show_header_program()
