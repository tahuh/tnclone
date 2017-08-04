#!/usr/bin/python
from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import SIGNAL
from PyQt4.Qt import *
import time, logging, StringIO
import os
import subprocess
from collections import Counter
import numpy as np
import operator
import sys

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
from Stats import ContigStats2 as ContigStats
from ErrorSeedCorrector import ErrorSeedCorrector
from Sampler import Sampler
from ReportError import ErrorReporter
from AlignmentChecker import SamReader , IndelMarker

CWD=os.getcwd()
"""
try:
	from Bio import SeqIO
	from Bio.Seq import Seq
	seqio_avail = True
except ImportError:
	from CustomParser import FastaParser
	from DNAUtils import translate
	seqio_avail = False
"""
from CustomParser import FastaParser
from SequenceParser import FAParser
from DNAUtils import translate

def abspath(path):
	return os.path.abspath(path)

def call_script(script):
	subprocess.call(script ,shell=True)
def launch(script):
	#commands.getoutput(script)
	subprocess.call(script, shell=True)
	
### Below Log file method from http://stackoverflow.com/questions/24469662/how-to-redirect-logger-output-into-pyqt-text-widget
# class QtHandler(logging.Handler):
	# def __init__(self):
		# logging.Handler.__init__(self)
	# def emit(self, record):
		# record = self.format(record)
		# if record : XStream.stdout().write("%s\n" %(record))
		
# logger = logging.getLogger(__name__)
# handler = QtHandler()
# handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
# logger.addHandler(handler)
# logger.setLevel(logging.DEBUG)

# class XStream(QtCore.QObject):
    # _stdout = None
    # _stderr = None
    # messageWritten = QtCore.pyqtSignal(str)
    # def flush( self ):
        # pass
    # def fileno( self ):
        # return -1
    # def write( self, msg ):
        # if ( not self.signalsBlocked() ):
            # self.messageWritten.emit(unicode(msg))
    # @staticmethod
    # def stdout():
        # if ( not XStream._stdout ):
            # XStream._stdout = XStream()
            # sys.stdout = XStream._stdout
        # return XStream._stdout
    # @staticmethod
    # def stderr():
        # if ( not XStream._stderr ):
            # XStream._stderr = XStream()
            # sys.stderr = XStream._stderr
        # return XStream._stderr
	
class logBuffer(QtCore.QObject, StringIO.StringIO):
    bufferMessage = QtCore.pyqtSignal(str)

    def __init__(self, *args, **kwargs):
        QtCore.QObject.__init__(self)
        StringIO.StringIO.__init__(self, *args, **kwargs)

    def write(self, message):
        if message:
            self.bufferMessage.emit(unicode(message))

        StringIO.StringIO.write(self, message)
		
		
class DiagnosisThreadPlot(QtCore.QThread):
	def __init__(self, data,canvas):
		QtCore.QThread.__init__(self)
		self.data = data
		self.canvas = canvas

	def plot(self):
		self.widget.canvas.ax.clear()
		self.widget.canvas.ax.plot(self.data)
		self.widget.canvas.ax.draw()


class DiagnosisThreadAnalyse(QtCore.QThread):
	def __init__(self,sample_info_file, trim_path, ref_path, editor):
		QtCore.QThread.__init__(self)
		self.sample_info = sample_info_file
		self.trim_path = trim_path
		self.ref_path = ref_path
		self.message = editor
		self.log_prefix = "[DiagnosisPlot] "
		self.root_path = "/".join(trim_path.split("/")[:-1])
	def __del__(self):
		self.wait()		
	def run(self):
		diagnosis_path = self.root_path + "/diagnosis"
		sam_diagnosis = diagnosis_path + "/sam"
		bam_diagnosis = diagnosis_path + "/bam"
		pileup_diagnosis = diagnosis_path + "/pileup"
		if not os.path.isdir(diagnosis_path):
			self.message.append(self.log_prefix + "Diagnosis folder do not exist. Create.")
			os.mkdir(diagnosis_path)

		if not os.path.isdir(sam_diagnosis):
			self.message.append(self.log_prefix + "Diagnosis->SAM folder do not exist. Create.")
			os.mkdir(sam_diagnosis)

		if not os.path.isdir(bam_diagnosis):
			self.message.append(self.log_prefix + "Diagnosis->BAM folder do not exist. Create.")
			os.mkdir(bam_diagnosis)

		if not os.path.isdir(pileup_diagnosis):
			self.message.append(self.log_prefix + "Diagnosis->Pileup folder do not exist. Create.")
			os.mkdir(pileup_diagnosis)
		cnt = 0
		with open(self.sample_info) as F:
			for line in F:
				data = line.rstrip().split()
				n = int(data[1])
				cnt += 1

		self.message.append(self.log_prefix + "Number of samples are {}. Diagnosis will take a while".format(cnt))
		t1 = time.time()
		with open(self.sample_info) as F :
			for line in F :
				data = line.rstrip().split("\t")
				gene = data[0]
				numbers = int(data[1])
				REF = self.ref_path + "/" + gene + ".gene.fa"
				bwa_index = "./bwa index %s"%(REF)
				self.message.append(self.log_prefix + bwa_index)
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
					bwa_mem = "./bwa mem -M -R '@RG\tID:%s\tSM:%s\tLB:%s\tPL:ILLUMINA\tPU:%s' %s %s %s > %s"%(gene,gene,gene,gene,REF,fwd_fastq,rev_fastq,bwa_sam)
					self.message.append(self.log_prefix + bwa_mem)
					call_script(bwa_mem)
					bam1 = bam_diagnosis + "/" + gene + "_" + str(i) + ".bam"
					samtools_view = "./samtools view -b -T %s -S %s -o %s"%(REF,bwa_sam,bam1)
					self.message.append(self.log_prefix + samtools_view)
					call_script(samtools_view)
					bam2 = bam1 + ".sorted"
					samtools_sort = "./samtools sort %s %s"%(bam1,bam2)
					self.message.append(self.log_prefix + samtools_sort)
					call_script(samtools_sort)
					#samtools_index = "./samtools index %s %s"%(bam2+ ".bam", bam2 + ".bam.bam")
					#self.message.append(self.log_prefix + samtools_index)
					#call_script(samtools_index)
					pileup_file = pileup_diagnosis + "/" + gene + "_" + str(i) + ".pileup"
					samtools_pileup = "./samtools mpileup -x -q 1 -Q 1 -d 1000000 -s -f %s %s > %s"%(REF, bam2+".bam",pileup_file)
					self.message.append(self.log_prefix + samtools_pileup)
					call_script(samtools_pileup)
		t2 = time.time()
		self.message.append(self.log_prefix + "Done diagnosis. Elapsed {} seconds".format(t2-t1))

### Main Thread

### Helper functions

def errorWindow(message):
   app = QApplication(sys.argv)
   w = QWidget()
   b = QPushButton(w)
   b.setText("TnClone Error Message")

   b.move(50,50)
   b.clicked.connect(lambda:showdialog(message))
   w.setWindowTitle("TnClone Error Window")
   w.show()
   sys.exit(app.exec_())
	
def showdialog(message):
   msg = QMessageBox()
   msg.setIcon(QMessageBox.Information)

   msg.setText("This is a message box")
   msg.setInformativeText("This is additional information")
   msg.setWindowTitle("MessageBox demo")
   msg.setDetailedText(message)
   msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)


class PipelineThread(QtCore.QThread):
	_log = QtCore.pyqtSignal(str)
	def __init__(self , opts,editor):
		QtCore.QThread.__init__(self)
		self.options = opts
		self.prefix = "TnClone"
		# self.prefix += time.strftime("%H-%M-%S")
		# self.prefix += "-"
		# self.prefix += time.strftime("%d-%m-%Y")
		self.message = editor
		self.trouble_shoot_file = open(self.options["out_path"] + "/troubleshoot.log" , "w")
		self.single_mode=False
		# XStream.stdout().messageWritten.connect(self.message.append)
		#self.log_ = QtCore.pyqtSignal(str)
	### Adding function required

	def __del__(self):
		self.wait()
	
	def flush_msg(self, module , record):
		#logger.info("[" + module + "] : " + record)
		_log_val = module +  " : " + record
		self._log.emit(_log_val)
		
		#c = self.message.textCursor()
		##print c
		#c.movePosition(QtGui.QTextCursor.End)
		#self.message.moveCursor(QtGui.QTextCursor.End)
		# self.display = QtGui.QTextBrowser()
		# self.display.verticalScrollBar().setValue(
		# self.display.verticalScrollBar().maximum())
		#self.message.ensureCursorVisible()
		#self.message.setTextCursor(c)
	def _run_sorter(self):
		outpath = abspath(self.options["out_path"])
		sort_path = outpath + "/sorted"
		logger = []
		#XStream.stdout().messageWritten.connect(self.message.append)
		if not os.path.isdir(sort_path) :
			#self.message.append("[LOGGER] Sorting file path does not exist. Creating....")
			#print ("TnClone::Sorter" , "Sorting file path does not exist. Creating....")
			#logging.info("[LOGGER] Sorting file does not exist. Creating....")
			os.mkdir(sort_path)

		if self.options["sort_on"] == True:
			sort_indir = os.path.abspath(self.options["raw_path"])
			sort_conf = abspath(self.options["sort_info"])
			sorter = Sorter(sort_indir,sort_path, sort_conf)
			# XStream.stdout().messageWritten.connect(self.message.append)
			#print ("TnClone::Sorter" , "Sorting samples...")
			#self.message.append("[LOGGER] Sorting your samples")
			#logging.info("[LOGGER] Sorting your samples")
#			try:
			sorter.sort()
			#print ("TnClone::Sorter" , "Done Sorting samples...")
#			except IOError:
#				errorWindow("TnClone has error in sorting. \nMay be there is no specified NGS file.\nOr check input file")

	def _run_trimmer(self):
		outpath = abspath(self.options["out_path"])
		sort_path = outpath + "/sorted"
		trim_path = outpath + "/trimmed"
		#XStream.stdout().messageWritten.connect(self.message.append)
		if not os.path.isdir(trim_path):
			#self.message.append("[LOGGER] Trimming file path does not exist. Creating....")
			#print ("TnClone::Trimmer" , "Trimming file path does not exist. Creating....")
			#logging.info("[LOGGER] Trimming file path does not exist. Creating....")
			os.mkdir(trim_path)

		trim_conf = abspath(self.options["sample_info"])

		trimmo_loc = abspath("./Trimmomatic-0.36/trimmomatic-0.36.jar")
		
		trimmer = Trimmer(sort_path, trim_path, trim_conf, trimmo_loc)
		#self.message.append("[LOGGER] Trimming your samples")
		#print ("TnClone::Trimmer" , "Trimming samples using Trimmomatic-0.36")
		#logging.info("[LOGGER] Trimming your samples")
		trimmer.trim()
		#print ("TnClone::Trimmer" , "Done trimming samples...")
	def _run_config_gen(self):
		
		logger = []

		outpath = self.options["out_path"]
		assem_conf_name = outpath + "/assembly.config"
		
		# if os.path.exists(assem_conf_name) : 
			# self.message.append("[LOGGER] Configuration file already exists. Program will not create.\n")
			# #logging.info("[LOGGER] Configuration file already exists. Program will not create.\n")
			# return
		# self.message.append("[LOGGER] Configuration file does not exist (AUTOMODE)\n")
		#XStream.stdout().messageWritten.connect(self.message.append)
		#self.message.append("[LOGGER] Building assembly configuration file")
		#print ("TnClone::ConfigGenerator" , "Building assembly configuration file")
		#logging.info("[LOGGER] Configuration file dows not exist (AUTOMODE)\n")
		#logging.info("[LOGGER] Configuration file dows not exist (AUTOMODE)\n")
		if self.options["from_bed"]:
			#self.message.append("[LOGGER] Generating from BED(region) file")
			for line in open(self.options["sample_info"]):
				data = line.rstrip().split()
				gene = data[0]
				nsamples = int(data[1])
				bed_path = self.options["bed_path"]
				ref_path = self.options["ref_path"]
				ksize = int(self.options["ksize"])
				cov = int(self.options["mink_occ"])
				##print gene , nsamples, bed_path, ref_path, ksize, cov
				script = "python buildConfigFile.py -k %d -b %s -c %d -n %d -r %s -o %s"%(ksize, bed_path, cov, nsamples, ref_path + "/" + gene + ".gene.fa", assem_conf_name)
				#self.message.append("[LOGGER] " + script)
				#print ("[TnClone::ConfigGenerator]" , script)
				#logging.info("[LOGGER] " + script)
				launch(script)

		else:
			#self.message.append("[LOGGER] Generating from specified source and terminal sequence")
			#print ("[TnClone::ConfigGenerator]" , "Not using region informaion. TnClone will directrly use sequence specified.")
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
		#XStream.stdout().messageWritten.connect(self.message.append)
		### Edit this

		### Required task

		### 1. Build Graph
		### 2. Modify graph
		### 2-1. If graph dissapears MARK as DISSAPEAR
		### 3. Assemble contigs
		### 3-1. If there is no seed match after error seed detection -> MARK as NO SEED and continue to next sample
		### 4. Select contigs with source and sink presents
		### If there is possible contigs 
		### 5. Pre-align
		### 6. If there is INDEL in contig, mark it as INDEL and keep store this sequence
		### Continue to assemble
		### If there is no possible contigs
		### MARK as BROKEN and continue
	
		### All the statistical analysis will be done in Report Analysis

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
				fwd_fqname = fq_path + "/" + gene + "/" + sample + "_R1.trimmed.fastq"
				rev_fqname = fq_path + "/" + gene + "/" + sample + "_R2.trimmed.fastq"

				fwd_size = os.stat(fwd_fqname).st_size
				rev_size = os.stat(rev_fqname).st_size

				down_sample = self.options["down_sample"]
				sampling_factor = float(self.options["down_sample_rate"])
				sampling = False
				if down_sample:
					if fwd_size / ( 1024 * 1024 ) >= 30 and rev_size / ( 1024 * 1024 ) >= 30 : ### Greater than 30Mb for uncompressed data
						sampling = True
						#print ("[TnClone::AssemblyMachinery]" , "Detected file size > 30mb... Downsampling process triggered...")
						#print ("[TnClone::AssemblyMachinery]" , "File size will be reduced from [FWD] %d Mb -> %d Mb || [REV] %d Mb -> %d Mb using factor %f\n"%(fwd_size / ( 1024 * 1024 ) , int(fwd_size / ( 1024 * 1024 ) * sampling_factor) , rev_size / ( 1024 * 1024 ) , int(rev_size / ( 1024 * 1024 ) * sampling_factor) , sampling_factor))
						#self.message.append("[LOGGER] Detected file size > 30mb... Downsampling process triggered...\n")
						#self.message.append("[LOGGER] File size will be reduced from [FWD] %d Mb -> %d Mb || [REV] %d Mb -> %d Mb using factor %f\n"%(fwd_size / ( 1024 * 1024 ) , int(fwd_size / ( 1024 * 1024 ) * sampling_factor) , rev_size / ( 1024 * 1024 ) , int(rev_size / ( 1024 * 1024 ) * sampling_factor) , sampling_factor))
						#logging.info("[LOGGER] Detected file size > 30mb... Downsampling process triggered...\n")
						#logging.info("[LOGGER] File size will be reduced from [FWD] %d Mb -> %d Mb || [REV] %d Mb -> %d Mb using factor %f\n"%(fwd_size / ( 1024 * 1024 ) , int(fwd_size / ( 1024 * 1024 ) * sampling_factor) , rev_size / ( 1024 * 1024 ) , int(rev_size / ( 1024 * 1024 ) * sampling_factor) , sampling_factor))
						sampler = Sampler(fwd_fqname,rev_fqname , sample_ratio=sampling_factor)
						sampler.sample()
					else:
						#self.message.append("[LOGGER] Eventhough down-sample is set, file size is smaller than 30 Mb. \nThis will not affect program performance thus there is no down sample process\n")
						#print ("[TnClone::AssemblyMachinery]" "Eventhough down-sample is set, file size is smaller than 30 Mb. \nThis will not affect program performance thus there is no down sample process")
						#logging.info("[LOGGER] Eventhough down-sample is set, file size is smaller than 30 Mb. \nThis will not affect program performance thus there is no down sample process\n")
						sampling = False
				if sampling:
					fwd_fqname = fq_path + "/" + gene + "/" + sample + "_R1.trimmed.fastq.dwn.fastq"
					rev_fqname = fq_path + "/" + gene + "/" + sample + "_R1.trimmed.fastq.dwn.fastq"

				#self.message.append("[LOGGER] Building KFS table....\n")
				#logging.info("[LOGGER] Building KFS table....\n")

				kfs_path = outpath + "/kfs"
				if not os.path.isdir(kfs_path):
					#self.message.append("[LOGGER] KFS does not exist. Create.")
					#print ("[TnClone::AssemblyMachinery]" , "KFS does not exist. Create.")
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

				# del srp
				# del kfs


				#self.message.append("[LOGGER] Build Graph using Illumina sequencing reads")
				#print ("[TnClone::AssemblyMachinery]" , "Build Graph using Illumina sequencing reads")
				tmp_path = outpath + "/tmp"

				if not os.path.isdir(tmp_path):
					#self.message.append("[LOGGER] temporal directory not exist. Create")
					#print ("[TnClone::AssemblyMachinery]" , "Temporal directory not exist. Create")
					os.mkdir(tmp_path)

				assem_path = outpath + "/assem"
				if not os.path.isdir(assem_path):
					#self.message.append("[LOGGER] No assembly result directory. Create...")
					#print ("[TnClone::AssemblyMachinery]" , "No assembly result directory. Create...")
					os.mkdir(assem_path)
				cfastq_path = outpath + "/cfastq"

				if not os.path.isdir(cfastq_path):
					#self.message.append("[LOGGER] No CFASTQ directory. Create")
					self.flush_smg("[TnClone::AssemblyMachinery]" , "No CFASTQ directory. Create")

					os.mkdir(cfastq_path)

				node_tmp = tmp_path + "/" + prefix + "_" + sample + "_k" + str(ksize) + ".nodes"

				assem_outname = assem_path + "/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig"
				mink_occ = int(self.options["mink_occ"])
				assembler = Assembler(kfs_oname, src, snk,  assem_outname, ksize, fwd_fqname, rev_fqname, node_tmp, mink_occ)

				#self.message.append("[LOGGER] Building....")
				#self.message.append("[MACHINERY] Building graph...\n")
				#print ("[TnClone::AssemblyMachinery]" , "Building graph....")
				
				assembler.buildGraph()
				
				#self.message.append("[LOGGER] Dump raw connection information....")
				#print ("[TnClone::AssemblyMachinery]" , "Dump connection information ( before modification )")

				assembler.dump2Tmp(node_tmp)
				#print ("[TnClone::AssemblyMachinery]" , "Modifying graph...")

				try:
					assembler.modifyGraph(node_tmp)
					graph_reduction_error = False
				except ValueError:
					#self.message.append("[LOGGER] Graph disappeared. TnClone cannot see ant graph that has connection info")
					#print ("[TnClone::AssemblyMachinery]" , "Graph disappeared. TnClone cannot see any graph that has connection info")
					self.trouble_shoot_file.write("At sample " + sample + ", graph information dissapeared. All output information for this sample will be marked as $")
					assem_out = open(assem_outname , "w")
					assem_out.write(">DISSAPEAR\nDISAPPEARED\n")
					assem_out.close()
					ss_find_file = open(outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.ss.rmv" , "w")
					ss_find_file.close()
					dup_rmv_file = open(outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.aln.target" , "w")
					dup_rmv_file.close()

					indel_tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.indel"
					good_tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.good"
					indel_f = open(indel_tmp_name , "w")
					good_f = open(good_tmp_name , "w")

					indel_f.close()
					good_f.close()
					continue



				#print ("[TnClone::AssemblyMachinery]" , "Dump connection information ( after modification )")
				assembler.dump2Tmp(node_tmp)
				
				
				#self.message.append("[LOGGER] Launching Assembly")
				#self.message.append("[MACHINERY] Assembly starts...\n")
				#print ("[TnClone::AssemblyMachinery]" , "Launching assembly core algorithm")
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
					#gpu_avail = self.options["gpu_avail"]

					ESC = ErrorSeedCorrector(G , seed, path, err_log, max_mis=max_mis,use_depth_profile=use_depth_correction)
					ESC.select()

				### PASTE HERE
				if assembly_seed_error:
					if os.stat(outpath + "/seed.fa").st_size == 0 :
						assem_out.write("NO SEED\n")
						assem_out.close()
						continue
					seed_file = open(outpath + "/seed.fa")
					chk = False
					for line in seed_file:
						src = line.rstrip()
						assembler.set_src(src)
						ret_val = self._assemble(assembler , assem_out, seed_error=True)
						if ret_val == 1 :
							break
				else:
					self._assemble(assembler, assem_out)

				assem_out.close()
				assembler.clearAssembler()

				#self.message.append("Finished assembly task for sample : %s"%(sample))
				
				### SEMI PROCESS
				ss_find_file = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.ss.rmv"
				ssf = SSF(assem_outname , ss_find_file, src, snk,self.message)
				ssf.run()

				dup_rmv_file = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.aln.target"
				drt = DRT(ss_find_file , dup_rmv_file,self.message)
				drt.run()

				### check if reference is available
				if self.options["ref_path"] == None or self.options["ref_path"] == '':
					refpath = False
				else:
					refpath = abspath(self.options["ref_path"])
				if refpath:
					### PRE-ALIGNMENT
					m = int(self.options["match_score"])
					x = int(self.options["mismatch_score"])
					g = int(self.options["gap_open"])
					e = int(self.options["gap_extension"])
					refname = refpath + "/" + sample.split("_")[0] + ".gene.fa"
					query = dup_rmv_file
					if not os.path.isdir(abspath(self.options["out_path"]) + "/pre-aln"):
						os.mkdir(abspath(self.options["out_path"]) + "/pre-aln")
					self.message.append("[TnClone::Pre-Alignment] Start...")
					samname = abspath(self.options["out_path"]) + "/pre-aln/" + sample + "_pre-aln.dna.sam"
					internal_align_script = "./SSW-CPP -m %d -x %d -g %d -e %d %s %s > %s"%(m,x,g,e,refname,query,samname)
					
					#self.message.append("[TnClone::Pre-Alignment] " + internal_align_script)
					#print ("[TnClone::Pre-Alignment]" , internal_align_script)
					
					launch(internal_align_script)
					#self.message.append("[TnClone::Pre-Alignment] Done pre-alignment...")

					### MARK UP INDELS
					sam_reader = SamReader(samname)
					sam_reader.open()
					#self.message.append("[TnClone::IndelMarker] Mark indel containing sequences")
					#print ("[TnClone::IndelMarker]" , "Mark indel containing sequences")
					marker = IndelMarker(sam_reader)
					marker.mark()

					#self.message.append("[TnClone::IndelMarker] Done marking")
					##print ("[TnClone::IndelMarker]" , "Done marking.")
					### TODO : IndelMark Flushing
					#self.message.append("[TnClone::AssemblyMachinery] Flushing...")					
					
					#print ("[TnClone::AssemblyMachinery]" , "Flushing marked sequences")
					cached = marker.cache
					ids = cached.id()
					seqs = cached.seq()
					indel_tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.indel"
					good_tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.good"
					indel_f = open(indel_tmp_name , "w")
					good_f = open(good_tmp_name , "w")
					
					Z = zip(ids, seqs)
					for id , seq in Z :
						if "INDEL" in id :
							indel_f.write(">" + id + "\n" + seq + "\n")
						else:
							good_f.write(">" + id + "\n" + seq + "\n")
					indel_f.close()
					good_f.close()
							
					#self.message.append("[TnClone::AssemblyMachinery] Done flushing...")
				else:
					### NO TASK ( MAFFT is too slow )
					indel_tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.indel"
					good_tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.good"
					indel_f = open(indel_tmp_name , "w")
					good_f = open(good_tmp_name , "w")
					indel_f.close()
					good_f.close()
					continue
		T2 = time.time()
		manager.profile()
		#self.message.append("[LOGGER] MEM used : " + manager.mem_usage + "\n")
		#self.message.append("[LOGGER] Performed in %f sec\n"%(T2-T1))
		
		#print ("[TnClone::AssemblyMachinery]", "MEM used for assembly : " + manager.mem_usage + "\n")
		#print ("[TnClone::AssemblyMachinery]", "Performed in %f sec\n"%(T2-T1))
	
	def _assemble(self, assembler,assem_out, seed_error=False):
		try:
			ncontig=0
			for path in assembler.assemble().walk():
				if "NOSEED" in path:
					if seed_error:
						return "NOSEED"
					assem_out.write("NO SEED\n")
					continue
				if "BROKEN" in path:
					niter = path[-1]
					total_depth , seq ,node_used, d_list , r_dict = path[0], path[1][0], path[2] , path[3], path[4]
					for entry in path[1][1:]:
						kmers = entry
						seq += kmers[-1]
						if seq != '':
							assem_out.write(">BROKEN|%d-th iteration\n"%(niter+1))
							assem_out.write(seq + "\n")

						else:				
							assem_out.write(">BROKEN|%d-th iteration\n"%(niter+1))
							assem_out.write("JUNKSEQUENCE\n")
							continue
				total_depth , seq ,node_used, d_list , r_dict = path[0], path[1][0], path[2] , path[3], path[4]
				r_list = r_dict.values()

				for entry in path[1][1:] :
					kmers = entry
					seq += kmers[-1]
				if seq != '' :
					assem_out.write(">contig:%d\n"%(ncontig))
					assem_out.write(seq + "\n")
					ncontig += 1
					continue

		except KeyError:
			
			### Internal Broken while traversing		
			assem_out.write(">BROKEN\n")
			assem_out.write("JUNKSEQUENCE\n")
				
		return 1

	def _run_alignment(self):
		#XStream.stdout().messageWritten.connect(self.message.append)
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
			#elf.message.append("[LOGGER] SAM folder do not exist. Create.")
			#print ("[TnClone::Alignment]" , "SAM folder do not exist. Create.")
			#logging.info("[LOGGER] SAM folder do not exist. Create.")
			os.mkdir(sam_path)
		outpath = abspath(self.options["out_path"])
		vcf_path = outpath + "/vcf"
		#log = []
		prefix = self.prefix
		if not os.path.isdir(vcf_path):
			#self.message.append("[LOGGER] No VCF dir. Create.")
			#print ("[TnClone::Aligment]" , "No VCF dir. Create.")
			#logging.info("[LOGGER] No VCF dir. Create.")
			os.mkdir(vcf_path)

		assem_fname = abspath(outpath + "/assembly.config")
		no_seed_samples = []
		broken_samples = []
		
		### These two containers will help user
		self.noseed_samples = []
		self.broken_samples = []
		self.graph_dissapeared_samples = []
		self.sample_counts_info = {}
		with open(assem_fname) as conf_file:
			for line in conf_file:
				data = line.rstrip().split("\t")
				sample = data[0]
				ksize = int(data[1])
				gene = sample.split("_")[0]
				if not self.sample_counts_info.has_key(gene):
					self.sample_counts_info[gene] = 1
				else:
					self.sample_counts_info[gene] += 1
				refname = refpath + "/" + sample.split("_")[0] + ".gene.fa"

				align_query = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.good"

				### MUST EDIT FROM HERE
				### Check if file size = 0

				fsize = os.stat(align_query).st_size
				if fsize == 0:
					### Figure out if the origin was no seed
					origin_file_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig"
					origin_file = open(origin_file_name)
					header = origin_file.readline()
					if header.startswith("NO SEED"):
						self.noseed_samples.append(sample)
						#continue
					elif header.startswith(">BROKEN"):
						self.broken_samples.append(sample)
						#continue
					elif header.startswith(">DISSAPEAR"):
						self.broken_samples.append(sample)
						#continue

				fsize2 = os.stat(outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.indel").st_size
				if fsize2 == 0 :
					tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.final.fa"
					tmp_merged = open(tmp_name , "w")
					tmp_merged.close()

				#self.message.append("[TnClone::Analysis] Merging indel and good files...\n")
				#print ("[TnClone::Analysis]" , "Merging indel and good files...")
				tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.final.fa"
				tmp_merged = open(tmp_name , "w")
				with open(align_query) as F:
					for line in F :
						tmp_merged.write(line)
				align_query = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.indel"
				with open(align_query) as F:
					for line in F:
						tmp_merged.write(line)
				tmp_merged.close()
				query = tmp_name

				### Perform alignment for query

				samname = sam_path + "/" + sample + "_k" + str(ksize) + ".dna.sam"

				script = "./SSW-CPP -m %d -x %d -g %d -e %d %s %s > %s"%(m,x,g,e,refname,query,samname)

				#self.message.append("[TnClone::Alignment] " + script+"\n")
				#print ("[TnClone::Alignment]" , script)
				launch(script)
				
				### enforce sleep for 0.5 sec to memory alignment
				QtCore.QThread.msleep(500)
				vcfname = vcf_path + "/" + sample + "_k" + str(ksize) + ".dna.vcf"
				script2 = "python variant_caller.py --sam %s --ref %s --vcf %s"%(samname, refname, vcfname)
				#self.message.append("[TnClone::VariantCaller] " + script2 + "\n")
				#print ("[TnClone::VariantCaller]", script2)
				
				launch(script2)
				aa_contigname = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.final.fa.aa.fa"
				aa_refname = refpath + "/" + prefix + "_" + sample.split("_")[0] + ".gene.fa.aa.fa"

				aa_contigf = open(aa_contigname , "w")
				aa_reff = open(aa_refname , "w")

				ref_parser = FAParser(refname)
				ref_parser.open()
				
				for id , desc , seq in ref_parser.parse():
					seq = seq.upper()
					aa = translate(seq)
					aa_reff.write(">"+id + "\n" + aa + "\n")
				aa_reff.close()
				ref_parser.close()
				# contig_parser = FastaParser(query)
				# contig_parser.open()
				# while 1:
					# try:
						# contig_parser.read()
						# cid = contig_parser.id
						# for frame in [0,1,2]:
							# caa = translate(contig_parser.seq[frame:])
							# aa_contigf.write(cid + "-frame:" + str(frame) + "\n" + caa + "\n")
					# except StopIteration:
						# break
				contig_parser = FAParser(query)
				contig_parser.open()
				for id , desc ,seq in contig_parser.parse():
					for frame in [0,1,2]:
						caa = translate(seq[frame:])
						aa_contigf.write(">" + id + "-frame:" + str(frame) + "\n" + caa + "\n")
						
				aa_contigf.close()
				contig_parser.close()

				acid_samname = sam_path + "/" + sample + "_k" + str(ksize) + ".protein.sam"
				script3 = "./SSW-C -m %d -x %d -o %d -e %d -p -a ./blosum62.txt -c -s -h %s %s > %s"%(m,x,g,e,aa_refname , aa_contigname, acid_samname)
				#self.message.append("[TnClone::Alignment] " + script3 + "\n")
				#print ("[TnClone::Alignment]" , script3)
				### enforce sleep for 0.5 sec to memory alignment
				QtCore.QThread.msleep(500)
				launch(script3)
				aa_vcfname = vcf_path + "/" + sample + "_k" + str(ksize) + ".protein.vcf"
				script4 = "python variant_caller.py --sam %s --ref %s --vcf %s"%(acid_samname, aa_refname, aa_vcfname)
				#self.message.append("[TnClone::VariantCaller] " + script4 + "\n")
				#print ("[TnClone::VariantCaller]" , script4)

				launch(script4)
	
	def _single_mode(self):
		m = int(self.options["match_score"])
		x = int(self.options["mismatch_score"])
		g = int(self.options["gap_open"])
		e = int(self.options["gap_extension"])
		try:
			outpath = abspath(self.options["out_path"])
			refpath = abspath(self.options["ref_path"])
		except KeyError:
			print "TnClone cannot configure mandatory fields."
			print "Please check if you pressed configure button"

		log = []

		sam_path = outpath + "/sam"
		prefix = self.prefix
		if not os.path.isdir(sam_path) :
			#elf.message.append("[LOGGER] SAM folder do not exist. Create.")
			#print ("[TnClone::Alignment]" , "SAM folder do not exist. Create.")
			#logging.info("[LOGGER] SAM folder do not exist. Create.")
			os.mkdir(sam_path)
		outpath = abspath(self.options["out_path"])
		vcf_path = outpath + "/vcf"
		#log = []
		prefix = self.prefix
		if not os.path.isdir(vcf_path):
			#self.message.append("[LOGGER] No VCF dir. Create.")
			#print ("[TnClone::Aligment]" , "No VCF dir. Create.")
			#logging.info("[LOGGER] No VCF dir. Create.")
			os.mkdir(vcf_path)

		assem_fname = abspath(outpath + "/assembly.config")
		no_seed_samples = []
		broken_samples = []
		
		### These two containers will help user
		self.noseed_samples = []
		self.broken_samples = []
		self.graph_dissapeared_samples = []
		self.sample_counts_info = {}
		
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
				fwd_fqname = fq_path + "/" + gene + "/" + sample + "_R1.trimmed.fastq"
				rev_fqname = fq_path + "/" + gene + "/" + sample + "_R2.trimmed.fastq"

				fwd_size = os.stat(fwd_fqname).st_size
				rev_size = os.stat(rev_fqname).st_size

				down_sample = self.options["down_sample"]
				sampling_factor = float(self.options["down_sample_rate"])
				sampling = False
				if down_sample:
					if fwd_size / ( 1024 * 1024 ) >= 30 and rev_size / ( 1024 * 1024 ) >= 30 : ### Greater than 30Mb for uncompressed data
						sampling = True
						#self.message.append("[LOGGER] Detected file size > 30mb... Downsampling process triggered...\n")
						#self.message.append("[LOGGER] File size will be reduced from [FWD] %d Mb -> %d Mb || [REV] %d Mb -> %d Mb using factor %f\n"%(fwd_size / ( 1024 * 1024 ) , int(fwd_size / ( 1024 * 1024 ) * sampling_factor) , rev_size / ( 1024 * 1024 ) , int(rev_size / ( 1024 * 1024 ) * sampling_factor) , sampling_factor))
						
						#self.flush_smg("[TnClone::AssemblyMachinery]" , "Detected file size > 30mb... Downsampling process triggered...\n")
						#print ("[TnClone::AssemblyMachinery]" , "File size will be reduced from [FWD] %d Mb -> %d Mb || [REV] %d Mb -> %d Mb using factor %f\n"%(fwd_size / ( 1024 * 1024 ) , int(fwd_size / ( 1024 * 1024 ) * sampling_factor) , rev_size / ( 1024 * 1024 ) , int(rev_size / ( 1024 * 1024 ) * sampling_factor) , sampling_factor))
						
						#logging.info("[LOGGER] Detected file size > 30mb... Downsampling process triggered...\n")
						#logging.info("[LOGGER] File size will be reduced from [FWD] %d Mb -> %d Mb || [REV] %d Mb -> %d Mb using factor %f\n"%(fwd_size / ( 1024 * 1024 ) , int(fwd_size / ( 1024 * 1024 ) * sampling_factor) , rev_size / ( 1024 * 1024 ) , int(rev_size / ( 1024 * 1024 ) * sampling_factor) , sampling_factor))
						sampler = Sampler(fwd_fqname,rev_fqname , sample_ratio=sampling_factor)
						sampler.sample()
					else:
						#self.message.append("[LOGGER] Eventhough down-sample is set, file size is smaller than 30 Mb. \nThis will not affect program performance thus there is no down sample process\n")
						#logging.info("[LOGGER] Eventhough down-sample is set, file size is smaller than 30 Mb. \nThis will not affect program performance thus there is no down sample process\n")
						sampling = False
				if sampling:
					fwd_fqname = fq_path + "/" + gene + "/" + sample + "_R1.trimmed.fastq.dwn.fastq"
					rev_fqname = fq_path + "/" + gene + "/" + sample + "_R1.trimmed.fastq.dwn.fastq"

				#self.message.append("[LOGGER] Building KFS table....\n")
				#print ("[TnClone::AssemblyMachinery]" , "Building KFS table....")
				#logging.info("[LOGGER] Building KFS table....\n")

				kfs_path = outpath + "/kfs"
				if not os.path.isdir(kfs_path):
					#self.message.append("[LOGGER] KFS does not exist. Create.")
					#print ("[TnClone::AssemblyMachinery]" , "KFS does not exist. Create.")
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

				# del srp
				# del kfs


				#print ("[TnClone::AssemblyMachinery]", "Build Graph using Illumina sequencing reads")

				tmp_path = outpath + "/tmp"

				if not os.path.isdir(tmp_path):
					#self.message.append("[LOGGER] temporal directory not exist. Create")

					os.mkdir(tmp_path)

				assem_path = outpath + "/assem"
				if not os.path.isdir(assem_path):
					#self.message.append("[LOGGER] No assembly result directory. Create...")

					os.mkdir(assem_path)
				cfastq_path = outpath + "/cfastq"

				if not os.path.isdir(cfastq_path):
					#self.message.append("[LOGGER] No CFASTQ directory. Create")

					os.mkdir(cfastq_path)

				node_tmp = tmp_path + "/" + prefix + "_" + sample + "_k" + str(ksize) + ".nodes"

				assem_outname = assem_path + "/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig"
				mink_occ = int(self.options["mink_occ"])
				assembler = Assembler(kfs_oname, src, snk,  assem_outname, ksize, fwd_fqname, rev_fqname, node_tmp, mink_occ)

				#self.message.append("[LOGGER] Building....")
				#print ("[TnClone::AssemblyMachinery]" , "Building graph..")

				
				assembler.buildGraph()
				
				#self.message.append("[LOGGER] Dump raw connection information....")

				assembler.dump2Tmp(node_tmp)
				#print ("[TnClone::AssemblyMachinery]", "Graph modification....")

				try:
					assembler.modifyGraph(node_tmp)
					graph_reduction_error = False
				except ValueError:
					#print ("[TnClone::AssemblyMachinery]", "[LOGGER] Graph disappeared. TnClone cannot see ant graph that has connection info")
					self.trouble_shoot_file.write("At sample " + sample + ", graph information dissapeared. All output information for this sample will be marked as $")
					assem_out = open(assem_outname , "w")
					assem_out.write(">DISSAPEAR\nDISAPPEARED\n")
					assem_out.close()
					ss_find_file = open(outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.ss.rmv" , "w")
					ss_find_file.close()
					dup_rmv_file = open(outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.aln.target" , "w")
					dup_rmv_file.close()

					indel_tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.indel"
					good_tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.good"
					indel_f = open(indel_tmp_name , "w")
					good_f = open(good_tmp_name , "w")

					indel_f.close()
					good_f.close()
					
					tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.final.fa"
					tmp_merged = open(tmp_name , "w")
					tmp_merged.close()
					
					continue



				#self.message.append("[LOGGER] Dump modified connection information....")
				assembler.dump2Tmp(node_tmp)
				
				
				#print ("[TnClone::AssemblyMachinery]", "Launching Assembly")
				#self.message.append("[MACHINERY] Assembly starts...\n")
				
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
					#gpu_avail = self.options["gpu_avail"]

					ESC = ErrorSeedCorrector(G , seed, path, err_log, max_mis=max_mis,use_depth_profile=use_depth_correction)
					ESC.select()

				### PASTE HERE
				if assembly_seed_error:
					if os.stat(outpath + "/seed.fa").st_size == 0 :
						assem_out.write("NO SEED\n")
						assem_out.close()
						
						ss_find_file = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.ss.rmv"
						ss = open(ss_find_file , "w")
						ss.close()
						
						dup_rmv_file = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.aln.target"
						dup = open(dup_rmv_file , "w")
						dup.close()
						
						indel_tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.indel"
						good_tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.good"
						indel_f = open(indel_tmp_name , "w")
						good_f = open(good_tmp_name , "w")

						indel_f.close()
						good_f.close()
						tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.final.fa"
						tmp_merged = open(tmp_name , "w")
						tmp_merged.close()
						continue
					seed_file = open(outpath + "/seed.fa")
					chk = False
					for line in seed_file:
						src = line.rstrip()
						assembler.set_src(src)
						ret_val = self._assemble(assembler , assem_out, seed_error=True)
						if ret_val == 1 :
							break
				else:
					self._assemble(assembler, assem_out)

				assem_out.close()
				assembler.clearAssembler()

				#self.message.append("Finished assembly task for sample : %s"%(sample))
				
				### SEMI PROCESS
				ss_find_file = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.ss.rmv"
				ssf = SSF(assem_outname , ss_find_file, src, snk,self.message)
				ssf.run()

				dup_rmv_file = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.aln.target"
				drt = DRT(ss_find_file , dup_rmv_file,self.message)
				drt.run()

				### check if reference is available
				if self.options["ref_path"] == None or self.options["ref_path"] == '':
					refpath = False
				else:
					refpath = abspath(self.options["ref_path"])
				if refpath:
					### PRE-ALIGNMENT
					m = int(self.options["match_score"])
					x = int(self.options["mismatch_score"])
					g = int(self.options["gap_open"])
					e = int(self.options["gap_extension"])
					refname = refpath + "/" + sample.split("_")[0] + ".gene.fa"
					#query = ss_find_file
					query = dup_rmv_file
					if os.stat(query).st_size == 0 :
						#print ("[TnClone::AssemblyMachinery]", "File size for " + dup_rmv_file + " is zero.")
						#print ("[TnClone::AssemblyMachinery]", "No pre-alignment sprocess.")
						#print ("[TnClone::AssemblyMachinery]", "Skip\n")
						indel_tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.indel"
						good_tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.good"
						indel_f = open(indel_tmp_name , "w")
						good_f = open(good_tmp_name , "w")

						indel_f.close()
						good_f.close()
						tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.final.fa"
						tmp_merged = open(tmp_name , "w")
						tmp_merged.close()
						
					if not os.path.isdir(abspath(self.options["out_path"]) + "/pre-aln"):
						os.mkdir(abspath(self.options["out_path"]) + "/pre-aln")
					#print ("[TnClone::Pre-Alignment]", "Start...")
					samname = abspath(self.options["out_path"]) + "/pre-aln/" + sample + "_pre-aln.dna.sam"
					if os.stat(query).st_size != 0 :
						
						internal_align_script = "./SSW-CPP -m %d -x %d -g %d -e %d %s %s > %s"%(m,x,g,e,refname,query,samname)
						# internal_align_script = "./SSW-C -m %d -x %d -o %d -e %d -c -s -h %s %s > %s"%(m,x,g,e,refname,query,samname)
						#print ("[TnClone::Pre-Alignment]",  internal_align_script)
						
						launch(internal_align_script)
						#self.message.append("[TnClone::Pre-Alignment] Done pre-alignment...")
						
						QtCore.QThread.msleep(5000)
						### MARK UP INDELS
						sam_reader = SamReader(samname)
						sam_reader.open()
						#print ("[TnClone::IndelMarker]","Mark indel containing sequences")
						marker = IndelMarker(sam_reader)
						marker.mark()

						#print ("[TnClone::IndelMarker]", "Done marking")

						### TODO : IndelMark Flushing
						#print ("[TnClone::AssemblyMachinery]" ,"Flushing...")					
						cached = marker.cache
						ids = cached.id()
						seqs = cached.seq()
						indel_tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.indel"
						good_tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.good"
						indel_f = open(indel_tmp_name , "w")
						good_f = open(good_tmp_name , "w")
						
						Z = zip(ids, seqs)
						for id , seq in Z :
							if "INDEL" in id :
								indel_f.write(">" + id + "\n" + seq + "\n")
							else:
								good_f.write(">" + id + "\n" + seq + "\n")
						indel_f.close()
						good_f.close()
						#self.message.append("[TnClone::AssemblyMachinery] Done flushing...")
					else:
						indel_tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.indel"
						good_tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.good"
						indel_f = open(indel_tmp_name , "w")
						good_f = open(good_tmp_name , "w")
						indel_f.close()
						good_f.close()
					
					#for line in conf_file:
					#	data = line.rstrip().split("\t")
					#	sample = data[0]
					#	ksize = int(data[1])
					#	gene = sample.split("_")[0]
					if not self.sample_counts_info.has_key(gene):
						self.sample_counts_info[gene] = 1
					else:
						self.sample_counts_info[gene] += 1
					refname = refpath + "/" + sample.split("_")[0] + ".gene.fa"

					align_query = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.good"

					
					### MUST EDIT FROM HERE
					### Check if file size = 0

					fsize = os.stat(align_query).st_size
					if fsize == 0:
						### Figure out if the origin was no seed
						origin_file_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig"
						origin_file = open(origin_file_name)
						header = origin_file.readline()
						
						if header.startswith("NO SEED"):
							self.noseed_samples.append(sample)
							#continue
						elif header.startswith(">BROKEN"):
							self.broken_samples.append(sample)
							#continue
						elif header.startswith(">DISSAPEAR"):
							self.broken_samples.append(sample)
							#continue
					fsize2 = os.stat(outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.indel").st_size
					if fsize2 == 0 :
						tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.final.fa"
						tmp_merged = open(tmp_name , "w")
						tmp_merged.close()
					
					#print ("[TnClone::Analysis]" , "Merging indel and good files...\n")
					tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.final.fa"
					tmp_merged = open(tmp_name , "w")
					with open(align_query) as F:
						for line in F :
							tmp_merged.write(line)
					align_query = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.indel"
					with open(align_query) as F:
						for line in F:
							tmp_merged.write(line)
					tmp_merged.close()
					query = tmp_name
					if os.stat(query).st_size != 0:
						### Perform alignment for query

						samname = sam_path + "/" + sample + "_k" + str(ksize) + ".dna.sam"

						script = "./SSW-CPP -m %d -x %d -g %d -e %d %s %s > %s"%(m,x,g,e,refname,query,samname)
						# script = "./SSW-C -m %d -x %d -o %d -e %d -c -s -h %s %s > %s"%(m,x,g,e,refname,query,samname)
						### enforce sleep for 5 sec to memory alignment
						#self.message.append("[TnClone::ToSystem] Wait for 5 sec to ensure memory to be aligned\n")
						QtCore.QThread.msleep(5000)
						#print ("[TnClone::Alignment]" , script+"\n")
						launch(script)
						vcfname = vcf_path + "/" + sample + "_k" + str(ksize) + ".dna.vcf"
						script2 = "python variant_caller.py --sam %s --ref %s --vcf %s"%(samname, refname, vcfname)
						#print ("[TnClone::VariantCaller]" , script2 + "\n")
						
						launch(script2)
						aa_contigname = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.final.fa.aa.fa"
						aa_refname = refpath + "/" + prefix + "_" + sample.split("_")[0] + ".gene.fa.aa.fa"

						aa_contigf = open(aa_contigname , "w")
						aa_reff = open(aa_refname , "w")
						
						# ref_parser = FastaParser(refname)
						# ref_parser.open()
						# ref_parser.read()
						# rid = ref_parser.id
						# aa = translate(ref_parser.seq)
						# aa_reff.write(rid + "\n" + aa + "\n")
						# aa_reff.close()
						# ref_parser.close()
						ref_parser = FAParser(refname)
						ref_parser.open()
						
						for id , desc , seq in ref_parser.parse():
							seq = seq.upper()
							aa = translate(seq)
							aa_reff.write(">" + id + "\n" + aa + "\n")
						aa_reff.close()
						ref_parser.close()
						# contig_parser = FastaParser(query)
						# contig_parser.open()
						# while 1:
							# try:
								# contig_parser.read()
								# cid = contig_parser.id
								# for frame in [0,1,2]:
									# caa = translate(contig_parser.seq[frame:])
									# aa_contigf.write(cid + "-frame:" + str(frame) + "\n" + caa + "\n")
							# except StopIteration:
								# break
						contig_parser = FAParser(query)
						contig_parser.open()
						for id , desc ,seq in contig_parser.parse():
							for frame in [0,1,2]:
								caa = translate(seq[frame:])
								aa_contigf.write(">"+id + "-frame:" + str(frame) + "\n" + caa + "\n")
								
						aa_contigf.close()
						contig_parser.close()

						acid_samname = sam_path + "/" + sample + "_k" + str(ksize) + ".protein.sam"
						script3 = "./SSW-C -m %d -x %d -o %d -e %d -p -a ./blosum62.txt -c -s -h %s %s > %s"%(m,x,g,e,aa_refname , aa_contigname, acid_samname)
						
						#self.message.append("[TnClone::Alignment] " + script3 + "\n")
						#print ("[TnClone::Alignment]" , script3 + "\n")

						launch(script3)
						### enforce sleep for 5 sec to memory alignment
						#self.message.append("[TnClone::ToSystem] Wait for 5 sec to ensure memory to be aligned\n")
						QtCore.QThread.msleep(5000)
						aa_vcfname = vcf_path + "/" + sample + "_k" + str(ksize) + ".protein.vcf"
						script4 = "python variant_caller.py --sam %s --ref %s --vcf %s"%(acid_samname, aa_refname, aa_vcfname)
						#print ("[TnClone::VariantCaller]" , script4 + "\n")

						launch(script4)
					else:
						continue ### No other task
				else:
					
					dup_rmv_file = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.aln.target"
					
					### Perform MAFFT
					
					cwd = os.getcwd()
					MAFFT_PATH = cwd + "/mafft"
					MAFFT_OUTPUT = outpath + "/sam/" + prefix + "_" + sample + "_k" + str(ksize) +".mafft"
					
					MAFFT_COMMAND = MAFFT_PATH + " " + dup_rmv_file  + " > " + MAFFT_OUTPUT
					
					launch(MAFFT_COMMAND)
					
					MSA2VCF_OUTPUT = outpath + "/vcf/" + prefix + "_" + sample + "_k" + str(ksize) +".dna.vcf"
					
					MSA2VCF_COMMAND = "java -jar " + cwd + "/msa2vcf.jar" + " " + MAFFT_OUTPUT + " > " + MSA2VCF_OUTPUT
					
					launch(MSA2VCF_COMMAND)
					
					tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.final.fa"
					tmp_merged = open(tmp_name , "w")
					with open(dup_rmv_file) as dup_tmp:
						for line in dup_tmp:
							tmp_merged.write(line)
					tmp_merged.close()
					
					indel_tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.indel"
					good_tmp_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.good"
					indel_f = open(indel_tmp_name , "w")
					good_f = open(good_tmp_name , "w")
					### We don't know real INDEL so we enforce all to analysis
					with open(dup_rmv_file) as dup_tmp:
						for line in dup_tmp:
							good_f.write(line)
					indel_f.close()
					good_f.close()
					
					continue
					
					
		
	def _perform_analysis(self):
		### Performing analysis using statistics called CCV or ratio


		def load_node(node_name):
			d = {}
			with open(node_name) as node_file:
				for line in node_file:
					data = line.rstrip().split("\t")
					node = data[0]
					depth = int(data[1])
					if len(data) == 2:
						d[node] = {"dp" : depth , "edge" :[] , "edge_dp" : [] }
					else:
						d[node] = {"dp" : depth , "edge" :[] , "edge_dp" : [] }
						edges = data[2].split(",")
						depths = map(lambda x : int(x) , data[3].split(","))
						d[node]["edge"] = edges
						d[node]["edge_dp"] = depths
			return d

		def load_snv_loc(contig_fname):
			### Contig_fname's contents must not contain INDEL containing sequence.
			var_pos = []
			small_parser = FAParser(contig_fname)
			small_parser.open()
			seqs = []
			for id , desc , seq in small_parser.parse():
				seqs.append(seq)
			Z = zip(*seqs)
			for i , z in enumerate(Z):
				c = Counter(z)
				if len(c) >1 :
					var_pos.append(i)

			return var_pos
		def load_sv_loc_denovo(vcf):
			var_pos = []
			with open(vcf) as VCF:
				for line in VCF:
					if line.startswith("#") : continue
					
					data = line.rstrip().split()
					reflen = len(data[3])
					altlen = len(data[4])
					
					if reflen >= 2 or altlen >= 2: continue ### Indel
					pos = int(data[1])-1
					var_pos.append(pos)
					
			return var_pos
		refpath = abspath(self.options["ref_path"])
		outpath = abspath(self.options["out_path"])
		assem_fname = abspath(outpath + "/assembly.config")
		self.ccv_pair = {}
		#self.message.append("[TnClone::Statistics] Computing Contig Coefficient of Variation(CCV) scores...")
		#print ("[TnClone::Statistics]" , "Computing Contig Coefficient of Variation(CCV) scores...")
		prefix = self.prefix
		tmp_path = outpath + "/tmp"
		with open(assem_fname) as conf_file:
			for line in conf_file:
				data = line.rstrip().split("\t")
				sample = data[0]
				ksize = int(data[1])
				target_query = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.good"
				node_name = tmp_path + "/" + prefix + "_" + sample + "_k" + str(ksize) + ".nodes"

				faparser = FAParser(target_query)
				faparser.open()
				cnt = 0
				for id , desc , seq in faparser.parse():
					cnt += 1
				faparser.close()
				if cnt == 1 :
					if os.stat(target_query).st_size == 0 :
						self.ccv_pair[sample] = {"n" : 0 , "ccv_set" : None}
						continue
					else:
						self.ccv_pair[sample] = {"n" : 1 , "ccv_set" : "s"}
						continue
				elif cnt == 0:
					self.ccv_pair[sample] = {"n" : 0 , "ccv_set" : None}
					continue
				else:
					### Possible analysis target
					### n describes how many contigs in SNV containing clone
					### 
					self.ccv_pair[sample] = {"n" : cnt , "ccv_set" : []}
					if refpath:
						snv_loc = load_snv_loc(target_query)
					else:
						vcfname = outpath + "/vcf/" + prefix + "_" + sample + "_k" + str(ksize) +".dna.vcf"
						snv_loc = load_snv_loc_denovo(vcfname)
					
					node_data = load_node(node_name)
					self.ccv_pair[sample]["ccv_set"] = self._compute_stat(target_query , node_data, snv_loc, cnt, ksize)
				
		#print ("[TnClone::Statistics]" ,"Done computing Contig Coefficient of Variation(CCV) scores...")
	def _get_confident_seqs(self):
		prefix = self.prefix
		#self.message.append("[TnClone::ContigSelector] Selecting Confident contigs using CCV value computed...")
		outpath = abspath(self.options["out_path"])
		assem_fname = abspath(outpath + "/assembly.config")
		with open(assem_fname) as conf_file:
			for line in conf_file:
				data = line.rstrip().split("\t")
				sample = data[0]
				ksize = int(data[1])

				target_query = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.good"
				indel_query = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.indel"
				final_file_name = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.final.confident.fa"
				final_file = open(final_file_name , "w")

				ccv_n = self.ccv_pair[sample]["n"]
				ccv_set = self.ccv_pair[sample]["ccv_set"]
				query_parser = FAParser()
				query_parser.open(target_query)

				if ccv_n == 0:
					# with open(indel_query) as F:
						# for line in F:
							# final_file.write(line)
					final_file.close()
					continue
				elif ccv_n == 1 :
					query_parser.close()
					with open(target_query) as F:
						for line in F:
							final_file.write(line)
					# with open(indel_query) as F:
						# for line in F:
							# final_file.write(line)
					final_file.close()
					continue
				elif ccv_n == 2 :
					ccv_ratio_map = dict(ccv_set)
					for id , desc , seq in query_parser.parse():
						id_new = id + "|ratio=" + str(ccv_ratio_map[id])
						final_file.write(">"+id_new + "\n" + seq + "\n")
					# with open(indel_query) as F:
						# for line in F:
							# final_file.write(line)
					final_file.close()
					continue
				else:
					ccv_ratio_map = dict(ccv_set)
					sorted_map = sorted(ccv_ratio_map.items() , key=operator.itemgetter(1))
					#selected = sorted_map[:2]
					#id_set = map(lambda x : x[0] , selected)
					selected = {}
					for id , value in sorted_map:
						### value is ccv
						if value <= 1.7:
							### Confident ones
							selected[id] = value
					id_set = map(lambda x : x[0] , selected)		
					for id , desc , seq in query_parser.parse():
						if not id in id_set : continue
						
						ccv = ccv_ratio_map[id]
						new_id = id + "|ccv=" + str(ccv)
						final_file.write(">"+new_id + "\n" + seq + "\n")
					# with open(indel_query) as F:
						# for line in F:
							# final_file.write(line)
					final_file.close()
					continue
					
		#self.message.append("[TnClone::ContigSelector] Done selection...")
	def _report_analysis(self):
		### Task
		### Load self.ccv_pair for count information (good one) -> got easier
		### Load final fa file to check how many indel seqs comparing n value above
		### Load SAM file to check if there is error in good contig DNA
		### Load SAM file to check if there is error in good contig protein


		authorship = """############################################
ANALYSIS INFORMATION
CREATED BY TnClone
PROGRAM AUTHOR : SUNGHOON HEO
PLEASE CONTACT AUTHOR VIA E-MAIL BELOW IF YOU HAVE ISSUE
team.tnclone(at)gmail.com
\n"""


		outpath = abspath(self.options["out_path"])
		assem_fname = abspath(outpath + "/assembly.config")
		prefix = self.prefix
		assem_path = outpath + "/assem"
		
		report_path = outpath + "/report"

		ksize = int(self.options["ksize"])


		if not os.path.isdir(report_path):
			#self.message.append("[LOGGER] No Report dir.Create.")
			#print ("[TnClone::Reporter]" , "No report directory. Create.")
			os.mkdir(report_path)

		report_name = outpath + "/report/analysis.summary"

		now = time.localtime()
		analysis_time = "%04d-%02d-%02d %02d:%02d:%02d" % (now.tm_year, now.tm_mon, now.tm_mday, now.tm_hour, now.tm_min, now.tm_sec)
		report = open(report_name , "w")
		report.write(authorship)
		report.write("ANALYSIS FINISHED AT : %s\n"%(analysis_time))
		report.write("###########################################\n")
#		report.close()

		#self.message.append("[TnClone::AnalysisReport] Reporting analysis result.\n")
		#print ("[TnClone::AnalysisReport]", "Reporting analysis result.\n")
#		confidence_report = open(report_path + "/dna_confident_contig_info.txt" , "w")
		count_table = {}
		for gene , n_sample in self.sample_counts_info.items():
			count_table[gene] = []
			for i in range(1,n_sample + 1):
				n_conf = 0
				tot_final_fa = 0
				sample = gene + "_" + str(i)
				indel_query = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.indel"
				parser = FAParser(indel_query)
				parser.open()
				n_indels = 0
				for id , desc , seq in parser.parse():
					n_indels += 1
				if os.stat(indel_query).st_size == 0 :
					n_indels = 0
				parser.close()
				
				final_fa = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.final.fa"
				if os.stat(final_fa).st_size == 0:
					tot_final_fa = 0
				else:
					tot_parser = FAParser(final_fa)
					tot_parser.open()
					for id , desc, seq in tot_parser.parse():
						tot_final_fa += 1
					tot_parser.close()
				confident_file = outpath + "/assem/" + prefix + "_" + sample + "_k" + str(ksize) + ".contig.final.confident.fa"
				#ccv_n = self.ccv_pair[sample]["n"]
				ccv_n = tot_final_fa
				if ccv_n >= 2:
					if os.stat(confident_file).st_size == 0 :
						n_conf = 0
					else:
						confident_parser = FAParser(confident_file)
						confident_parser.open()
						for id , desc , seq in confident_parser.parse():
							n_conf += 1
							
						confident_parser.close()
					total_count = "%d(C=%d:NC=%d:I/D=%d)"%(tot_final_fa , n_conf , tot_final_fa - n_conf, n_indels)
				else:
					if os.stat(confident_file).st_size == 0 :
						n_conf = 0
					else:
						confident_parser = FAParser(confident_file)
						confident_parser.open()
						for id , desc , seq in confident_parser.parse():
							n_conf += 1
							
						confident_parser.close()
					
					total_count = "%d(C=%d:I/D=%d)"%(tot_final_fa,n_conf, n_indels)
				
				count_table[gene].append(total_count)
#				confidence_report.write(sample)
#				ccv_set = self.ccv_pair[sample]
				
		#confidence_report.close()
		report.write("[Contig Number Summary]\n")
		for gene in count_table:
			report.write(gene + "\t")
			tabbed_str = "\t".join(count_table[gene])
			report.write(tabbed_str + "\n")

		### Loading DNA sam file
		dna_error_table = {}
		sam_path = outpath + "/sam"

		dna_free_report = open(report_path + "/dna_free_contig_info.txt" , "w")
		
		free_query = {}
		for gene , n_samples in self.sample_counts_info.items():
			dna_error_table[gene] = []
			free_query[gene] = {}
			for i in range(1 , n_samples + 1):
				sample = gene + "_" + str(i)
				samname = sam_path + "/" + sample + "_k" + str(ksize) + ".dna.sam"
				free_query[gene][i] = []
				if os.path.exists(samname) :
					sam_reader = SamReader(samname)
					sam_reader.open()
					form = []
					
					for record in sam_reader.read_record():
						cigar = record.cigar
						alphas = sam_reader.cigar_chars(cigar)
						query = record.query
						if len(alphas) == 1:
							if alphas[0] == 'M' or alphas[0] =='=' :
								free_query[gene][i].append(query)
								form.append("T")
							else:
								form.append("E")
						else:
							form.append("E") 
					sam_reader.close()
					converted_form = "(" + "/".join(form) + ")"
					dna_error_table[gene].append(converted_form)
				else:
					dna_error_table[gene].append("(NOINFO)")

		report.write("\n[DNA ERROR/FREE SUMMARY]\n")
		for gene in dna_error_table:
			report.write(gene + "\t")
			st_form = "\t".join(dna_error_table[gene])
			report.write(st_form + "\n")
			N = self.sample_counts_info[gene]
			dna_free_report.write(gene +"\t")
			for i in range(1,N+1):
				free_contigs = free_query[gene][i]
				if len(free_contigs) == 0:
					dna_free_report.write("\tNONE")
				else:
					dna_free_report.write( "\t" + ",".join(free_contigs) )
			dna_free_report.write("\n")


		### Loading protein SAM file
		aa_free_report = open(report_path + "/protein_free_contig_info.txt" , "w")
		aa_error_table = {}
		free_query = {}
		for gene , n_samples in self.sample_counts_info.items():
			aa_error_table[gene] = []
			free_query[gene] = {}
			for i in range(1 , n_samples + 1):
				sample = gene + "_" + str(i)
				samname = sam_path + "/" + sample + "_k" + str(ksize) + ".protein.sam"
				free_query[gene][i] = []
				if os.path.exists(samname):
				
					sam_reader = SamReader(samname)
					sam_reader.open()
					form = []

					for record in sam_reader.read_record():
						cigar = record.cigar
						alphas = sam_reader.cigar_chars(cigar)
						query = record.query
						if len(alphas) == 1:
							if alphas[0] == 'M' or alphas[0] =='=' :
								free_query[gene][i].append(query)
								form.append("T")
							else:
								form.append("E")
						else:
							form.append("E") 
					sam_reader.close()
					converted_form = "(" + "/".join(form) + ")"
					aa_error_table[gene].append(converted_form)
				else:
					aa_error_table[gene].append("(NOINFO)")

		report.write("\n[PROTEIN ERROR/FREE SUMMARY]\n")
		for gene in aa_error_table:
			report.write(gene + "\t")
			st_form = "\t".join(aa_error_table[gene])
			report.write(st_form + "\n")
			N = self.sample_counts_info[gene]
			aa_free_report.write(gene +"\t")
			for i in range(1,N+1):
				free_contigs = free_query[gene][i]
				if len(free_contigs) == 0:
					aa_free_report.write("\tNONE")
				else:
					aa_free_report.write( "\t" + ",".join(free_contigs) )
			aa_free_report.write("\n")
		#self.message.append("[TnClone::AnalysisReport] Reporting analysis done.")
		#print ("[TnClone::AnalysisReport]" , "Reporting analysis done.")
		#self.message.append("\n#####################################\n")
		#print ("[TnClone::AnalysisReport]" , "#####################################")
		#self.message.append("[TnClone] Finished whole analysis. Please look at results at " + outpath)
		#print ("[TnClone]" , "Finished whole analysis. Please look at results at " + outpath)
		#print ("[TnClone::AnalysisReport]" , "\n#####################################\n")
		#self.message.append("\n#####################################\n")

		report.close()
		dna_free_report.close()
		aa_free_report.close()

	def _compute_stat(self, contig_fname , node_data , var_pos, cnt , ksize):
		parser = FAParser(contig_fname)
		parser.open()
		ret = []
		for id ,desc , seq in parser.parse():
			ratio_list = []
			for p in var_pos:
				if p < ksize:
					kmer_before = seq[p-ksize:p]
					next_alleles = node_data[kmer_before]["edge"]

					other_kmers = map(lambda x : kmer_before[1:] + x , next_alleles)

					my_kmer = kmer_before[1:] + seq[p]

					my_depth = node_data[my_kmer]["dp"]

					other_depths = map(lambda x : node_data[x]["dp"] , other_kmers)

					norm_ratio = my_depths / float(sum(other_depths))
					ratio_list.append(norm_ratio)
			
			if len(ratio_list) == 0:
				ret.append( ( id, 'nan') )
			else:
				mean = np.mean(ratio_list)
				std = np.std(ratio_list)
				cv = std / mean
				if cnt > 2 :
					ret.append( ( id , cv ) )
				else:
					ret.append( ( id , mean ) )
		parser.close()
		return ret
	
	def run(self):
		#if not self.is_running:
		whole_start = time.time()
		if self.options["sort_on"] == True:
			print "[TnClone::Sorter] Demultiplexing samples...."
			self._run_sorter()
		if self.options["trim_on"] == True:
			print "[TnClone::Trimmer] Trimming samples...."
			self._run_trimmer()
		
		if self.options["assembly_on"] == True:
			print "[TnClone::ConfigGenerator] Generating assembly configuration"
			self._run_config_gen()
			print "[TnClone::AssemblyMachinery] Machinery starts... Using single mode."
			#self._run_until_assembly()
			### Single clone computing mode
		
			self.single_mode = True
			self._single_mode()
		if self.options["analysis_on"] == True:
			if self.options["ref_path"]== None or self.options["ref_path"] == '':
				pass
			else:
				if self.single_mode == True:
					print "[TnClone::DownstreamAnalysis] Performing analysis of assembled contigs...."
					self._perform_analysis()
					print "[TnClone::DownstreamAnalysis] Get confident seuqneces from assembled contigs...."
					self._get_confident_seqs()
					print "[TnClone::DownstreamAnalysis] Reporting analysis of assembled contigs...."
					self._report_analysis()
					print "[TnClone::DownstreamAnalysis] Done downstream analysis of assembled contigs...."
					print "[TnClone::DownstreamAnalysis] Result can be found at " , self.options["out_path"]
				else:
					print "[TnClone::DownstreamAnalysis] Start alignment of assembled contigs..."
					self._run_alignment()
					print "[TnClone::DownstreamAnalysis] Performing analysis of assembled contigs...."
					self._perform_analysis()
					print "[TnClone::DownstreamAnalysis] Get confident seuqneces from assembled contigs...."
					self._get_confident_seqs()
					print "[TnClone::DownstreamAnalysis] Reporting analysis of assembled contigs...."
					self._report_analysis()
					print "[TnClone::DownstreamAnalysis] Done downstream analysis of assembled contigs...."
					print "[TnClone::DownstreamAnalysis] Result can be found at " , self.options["out_path"]
		whole_end = time.time()
		print "\n===============================================================\n"
		print "[TnClone] Finished analysis... One can quit this appication."
		print "[TnClone] Elapsed " , (whole_end - whole_start) , " seconds for given analysis steps."
		print "\n===============================================================\n"
		#self.message.append("\n\nELAPSED TIME %f sec\n"%(whole_end - whole_start))
			#logging.info("\n\nELAPSED TIME %f sec\n"%(whole_end - whole_start))
		#else:
			#self.run()
	# def stop(self):
		# self.is_running = False
		# self.quit()
