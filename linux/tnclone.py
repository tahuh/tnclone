#!/usr/bin/python

import sys
import time
sys.stdout.write("sys module import\n")
from PyQt4 import QtGui, QtCore
sys.stdout.write("Qt -1\n")
from PyQt4.QtCore import SIGNAL
sys.stdout.write("Qy -2\n")
from PyQt4.Qt import *
sys.stdout.write('Qt -3\n')
from Queue import Queue
sys.stdout.write('Queue\n')
import window
sys.stdout.write('wndow\n')
#print "loaded custom wnidows"

import os
import commands
import subprocess

try:
	from Bio import SeqIO
	from Bio.Seq import Seq
	seqio_avail = True
except ImportError:
	from CustomParser import FastaParser
	from DNAUtils import translate
	seqio_avail = False

import logging


### CUSTOM module (SCREENS)
from PipelineThread import *
from PopupScreens import *



### MISC

def abspath(path):
	return os.path.abspath(path)

def launch(script):
	#commands.getoutput(script)
	subprocess.call(script, shell=True)

### Option class
class Option(object):
	def __init__(self, option_dict=None):
		self.option = option_dict





### Main UI class
class MainClass(QtGui.QMainWindow, window.Ui_MainWindow):
	def __init__(self):
		super(self.__class__,self).__init__()
		self.setupUi(self)

		self.user_default_dir = os.path.expanduser('~')
		### Required option fields
		self.pushButton_5.clicked.connect(self.browse_sortinfo)
		self.pushButton_6.clicked.connect(self.browse_sampleinfo)
		self.pushButton_4.clicked.connect(self.browse_rawpath)
		self.pushButton_10.clicked.connect(self.browse_refpath)
		self.pushButton_7.clicked.connect(self.browse_outpath)

		### Popup Options
		self.pushButton_9.clicked.connect(self.showPopupScreen)
		self.pushButton_11.clicked.connect(lambda:self.show_format(desc_type="sort"))
		self.pushButton_12.clicked.connect(lambda:self.show_format(desc_type="sample"))
		self.pushButton_13.clicked.connect(lambda:self.show_format(desc_type="fasta"))

		### Command buttons
		self.pushButton_2.clicked.connect(self.quit_application)
		self.pushButton.clicked.connect(self.abort_application)
		self.pushButton_3.clicked.connect(self.run_application)
		self.pushButton_8.clicked.connect(self.configure_application)

		self.pushButton_14.clicked.connect(self.showDiagnosisPlot)
		"""
		self.logBuffer = logBuffer()
		self.logBuffer.bufferMessage.connect(self.run_application)
		logFormatter = logging.Formatter('%(levelname)s: %(message)s')
		logHandler = logging.StreamHandler(self.logBuffer)
		logHandler.setFormatter(logFormatter)
		self.logger = logging.getLogger()
		self.logger.setLevel(logging.INFO)
		self.logger.addHandler(logHandler)
		"""
		
		self.popup_screen = None
		self.diagnosis_plot = None
		self.helper = None
		self.option_dict = {}
		self.kmer_update = '63'
		self.mink_occ = '3'
		self.sort_on = True
		self.trim_on = True
		self.analysis_on = True
		self.assembly_on = True
		self.bed_path = None
		self.match_score = '1'
		self.mismatch_score = '3'
		self.gap_open = '5'
		self.gap_extension = '3'
		self.down_sample = True
		self.down_sample_rate = '0.3'
		self.depth_correction = False
		self.mismatch_correction = False
		self.max_mis = '3'
		self.gpu_avail = False
		self.log_screen = self.textEdit

		self.raw_path = None
		self.sort_info = None
		self.sample_info = None
		self.ref_path = None
		self.source_seq = None
		self.sink_seq = None
		self.output_path = None
		self.config_from_bed = False
		self.bed_path = None

		### Some assertion field if any one of enry here false, program will not run
		self.kmer_length_assert = False
		self.bed_assert = False
		self.sort_info_assert = False
		self.sample_info_assert = False
		self.ref_path_assert = False
		self.raw_path_assert = False
		self.output_path_assert = False
		self.file_path = None
	### Show popup screen window
	def showPopupScreen(self):
		self.popup_screen = OptionPopupScreen()
		self.popup_screen.kmer_update.connect(self.changeProgValKmer)
		self.popup_screen.mink_occ.connect(self.changeProgValMinOcc)
		self.popup_screen.sort_on.connect(self.changeProgValSort)
		self.popup_screen.trim_on.connect(self.changeProgValTrim)
		self.popup_screen.analysis_on.connect(self.changeProgValAnal)
		self.popup_screen.assembly_on.connect(self.changeProgValAssem)
		self.popup_screen.bed_path.connect(self.changeProgValBed)
		self.popup_screen.match_score.connect(self.changeProgValMatch)
		self.popup_screen.mismatch_score.connect(self.changeProgValMisMatch)
		self.popup_screen.gap_open.connect(self.changeProgValGapO)
		self.popup_screen.gap_extension.connect(self.changeProgValGapE)

		self.popup_screen.down_sample.connect(self.changeProgValDownSample)
		self.popup_screen.down_sample_rate.connect(self.changeProgValDownSampleRate)
		self.popup_screen.depth_correction.connect(self.changeProgValDepthCorr)
		self.popup_screen.mismatch_correction.connect(self.changeProgValMismatch)
		self.popup_screen.max_mis.connect(self.changeProgValMaxMis)
#		self.popup_screen.gpu_avail.connect(self.changeProgValGpu)

		self.popup_screen.show()
		
		self.program = None
	def showDiagnosisPlot(self):
		if self.sample_info == None:
			self.log_screen.append("[LOGGER] Please trim your sample and insert sample information field first....")
			if self.ref_path == None:
				self.log_screen.append("[LOGGER] Please reference field first....")
			pass

		trim_path = self.output_path + "/trimmed"
		ref_path = self.ref_path 
		self.diagnosis_plot = DiagnosisPlot(self.sample_info,trim_path, ref_path)
		self.diagnosis_plot.show()

	def changeProgValDownSample(self,value):
		self.down_sample = value
	def changeProgValDownSampleRate(self,value):
		self.down_sample_rate = value
	def changeProgValDepthCorr(self, value):
		self.depth_correction = value
	def changeProgValMismatch(self, val):
		self.mismatch_correction = val
	def changeProgValMaxMis(self, value):
		self.max_mis = value
	def changeProgValGpu(self, val):
		self.gpu_avail = val
	#### Update option changes
	def changeProgValKmer(self,value):
		self.kmer_update = value
	def changeProgValMinOcc(self,value):
		self.mink_occ = value
	def changeProgValSort(self,value):
		self.sort_on = value
	def changeProgValTrim(self,value):
		self.trim_on = value
	def changeProgValAnal(self,value):
		self.analysis_on = value
	def changeProgValAssem(self,value):
		self.assembly_on = value
	def changeProgValBed(self,value):
		self.bed_path = value
	def changeProgValMatch(self, value):
		self.match_score = value
	def changeProgValMisMatch(self,value):
		self.mismatch_score = value
	def changeProgValGapO(self, value):
		self.gap_open = value
	def changeProgValGapE(self, value):
		self.gap_extension = value

	### Browsings 
	def browse_rawpath(self) :
		if self.file_path == None :
			filepath = QtGui.QFileDialog.getExistingDirectory(self.centralwidget, "Browse RAW NGS sequencing path",self.user_default_dir, QtGui.QFileDialog.ShowDirsOnly)
			self.file_path = str(filepath)
		else:
			filepath = QtGui.QFileDialog.getExistingDirectory(self.centralwidget, "Browse RAW NGS sequencing path",self.file_path, QtGui.QFileDialog.ShowDirsOnly)
		
		self.raw_path = str(filepath)
		self.lineEdit_7.setText(filepath)

	def browse_sortinfo(self) :
		if self.file_path == None:

			sort_info_path = QtGui.QFileDialog.getOpenFileName(self.centralwidget, "Browse Sorting config file" , self.user_default_dir , "*")
			self.file_path = str(sort_info_path)
		else:
			sort_info_path = QtGui.QFileDialog.getOpenFileName(self.centralwidget, "Browse Sorting config file" , self.file_path , "*")
		self.sort_info = str(sort_info_path)
		self.lineEdit_8.setText(sort_info_path)

	def browse_sampleinfo(self) :
		if self.file_path == None:
			sample_info_path = QtGui.QFileDialog.getOpenFileName(self.centralwidget,"Browse Sample info file" , self.user_default_dir , "*")
			self.file_path = sample_info_path
		else:
			sample_info_path = QtGui.QFileDialog.getOpenFileName(self.centralwidget,"Browse Sample info file" , self.file_path , "*")
			
		self.sample_info = str(sample_info_path)
		self.lineEdit_9.setText(sample_info_path)

	def browse_refpath(self):
		if self.file_path == None:
			refpath = QtGui.QFileDialog.getExistingDirectory(self.centralwidget , "Browse Reference location" , self.user_default_dir,QtGui.QFileDialog.ShowDirsOnly)
			self.file_path = str(refpath)
		else:
			refpath = QtGui.QFileDialog.getExistingDirectory(self.centralwidget , "Browse Reference location" , self.file_path,QtGui.QFileDialog.ShowDirsOnly)
		self.ref_path = str(refpath)
		self.lineEdit_14.setText(refpath)

	def browse_outpath(self):
		if self.file_path == None:
			filepath = QtGui.QFileDialog.getExistingDirectory(self.centralwidget, "Browse output directory",self.user_default_dir,QtGui.QFileDialog.ShowDirsOnly)
			self.file_path = str(file_path)
		else:
			filepath = QtGui.QFileDialog.getExistingDirectory(self.centralwidget, "Browse output directory",self.file_path,QtGui.QFileDialog.ShowDirsOnly)

		self.output_path = str(filepath)
		self.lineEdit_12.setText(filepath)
	### Format helper showing script
	def show_format(self, desc_type=None):
		self.helper = DescriptionPopupScreen(desc_type=desc_type)
		self.helper.setGeometry(QRect(0,0,1000,600))
		self.helper.setWindowTitle(desc_type.upper() + " Information")
		self.helper.show()

	def configure_application(self):

		### Default options on main screen

		self.option_dict["ksize"] = int(self.kmer_update)
		self.option_dict["mink_occ"] = int(self.mink_occ)

		self.option_dict["raw_path"] = self.raw_path

		if self.raw_path == None:
			self.raw_path_assert = None

		self.option_dict["sort_info"] = self.sort_info

		if self.sort_info == None or self.sort_info == '' and self.sort_on == True:
			self.sort_info_assert = True

		self.option_dict["sample_info"] = self.sample_info

		if self.sample_info == None or self.sample_info == '' and self.trim_on == True:
			self.sample_info_assert = True

		self.source_seq = str(self.lineEdit_10.text())
		self.sink_seq = str(self.lineEdit_11.text())
		self.option_dict["source_seq"] = self.source_seq
		self.option_dict["sink_seq"] = self.sink_seq
		
		if self.source_seq == None or self.sink_seq == None :
			self.config_from_bed = True
			
		### Check if bed file avail if not then Tell log message
		self.option_dict["from_bed"] = False
		if self.bed_path == None and self.config_from_bed == True:
			self.bed_assert = True
		if self.bed_path != '' and self.bed_path != None:
			self.option_dict["from_bed"] = True

		if self.config_from_bed == False and self.bed_path == None:

			if len(self.source_seq) != int(self.kmer_update):
				self.kmer_length_assert = True
			if len(self.sink_seq) != int(self.kmer_update):
				self.kmer_length_assert = True
		self.option_dict["bed_path"] = self.bed_path
		self.option_dict["ref_path"] = self.ref_path
		if self.ref_path == None:
			self.ref_path_assert = None



		self.option_dict["out_path"] = self.output_path
		if self.output_path == None:
			self.output_path_assert == True

		self.option_dict["sort_on"] = self.sort_on
		self.option_dict["trim_on"] = self.trim_on
		self.option_dict["assembly_on"] = self.assembly_on
		self.option_dict["analysis_on"] = self.analysis_on

		### Options in pop up screen
		self.option_dict["match_score"] = self.match_score
		self.option_dict["mismatch_score"] = self.mismatch_score
		self.option_dict["gap_open"] = self.gap_open
		self.option_dict["gap_extension"] = self.gap_extension

		"""if self.kmer_length_assert == False and \
			self.bed_assert == False and \
			self.sort_info_assert == False and \
			self.sample_info_assert == False and \
			self.ref_path_assert == False and \
			self.raw_path_assert == False and \
			self.output_path_assert == False :
			self.show_configuration()
		"""
		#if self.raw_path_assert == False:
		#	self.show_configuration()
		self.show_configuration()

		self.option_dict["down_sample"] = self.down_sample
		self.option_dict["down_sample_rate"] = self.down_sample_rate
		self.option_dict["depth_correction"] = self.depth_correction
		self.option_dict["mismatch_correction"] = self.mismatch_correction
		self.option_dict["max_mis"] = self.max_mis
		#self.option_dict["gpu_avail"] = self.gpu_avail
	def show_configuration(self):
		string = """\n[CONFIGURATION INFORMATION]
############################################################
ASSEMBLY OPTIONS

NGS path              : {}
K-mer size            : {}
Min k-mer occ         : {}
Start sequence        : {}
End sequence          : {}
Sort information      : {}
Sample information    : {}

ALIGNMENT OPTIONS

Reference directory   : {}
BED file location     : {}
Match score           : {}
Mismatch score        : {}
Gap open penalty      : {}
Gap extension penalty : {}

OTHERS

Output directory      : {}
Analysis Steps 
         Sort         : {}
         Trim         : {}
         Assembly     : {}
         DownStream   : {}

Down Sample           : {}
Down sample rate      : {}
Depth correction      : {}
Mismatch correction   : {}
Number of mismatched  : {}
############################################################\n
""".format(self.raw_path,self.kmer_update, self.mink_occ,  self.source_seq, self.sink_seq, self.sort_info, self.sample_info,self.ref_path,self.bed_path,self.match_score, self.mismatch_score, self.gap_open, self.gap_extension, self.output_path,self.sort_on, self.trim_on, self.assembly_on, self.analysis_on,self.down_sample, self.down_sample_rate, self.depth_correction, self.mismatch_correction, self.max_mis)
		self.log_screen.append(string)
	### Misc
	def set_source(self) :
		self.source_seq = str(self.lineEdit_10.text())
		if self.source_seq == None or self.source_seq == '':
			self.source_seq = None
	def set_sink(self):
		self.sink_seq = str(self.lineEdit_11.text())
		if self.sink_seq == None or self.sink_seq == '':
			self.sink_seq = None
	@pyqtSlot()
	def run_application(self):
		self.program = PipelineThread(self.option_dict, self.textEdit)
		
		#self.program.start()
		self.thread = QtCore.QThread()
		self.program = PipelineThread(self.option_dict, self.textEdit)
		self.program.moveToThread(self.thread)
		self.thread.started.connect(self.program.run)
		self.thread.start()
		self.connect(self.program,SIGNAL('finished()'),self.done)
		#self.connect(self.thread,SIGNAL(('finished()'),self.done))
	@pyqtSlot(str)
	def append_text(self, text):
		self.textEdit.moveCursor(QtGui.QTextCursor.End)
		self.textEdit.insertPlainText(text)

	def done(self):
		self.program.quit()
		QtGui.QMessageBox.information(self, "Done!" , "Done Analysis!!")
		self.program.wait()

	def abort_application(self):
		choice = QtGui.QMessageBox.question(self.centralwidget, "Abort program", "Are you sure to abort?" , QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)

		if choice == QtGui.QMessageBox.Yes:
#			self.pipe_thread.stop()
			self.program.quit()
			self.program.wait()
			QtGui.QMessageBox.information(self, "Abort!", "Aborted...")
			#print "Aborting..."
		else:
			pass
	def quit_application(self):
		choice = QtGui.QMessageBox.question(self.centralwidget,"Exit Program", "Are you sure to exit the program?", QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
		if choice == QtGui.QMessageBox.Yes:
			sys.exit()
		else:
			pass

class WriteStream(object):
	def __init__(self, queue):
		self.queue = queue
	def write(self, text):
		self.queue.put(text)

class MyReciever(QObject):
	mySignal = QtCore.pyqtSignal(str)
	def __init__(self,queue,*args,**kwargs):
		QObject.__init__(self,*args,**kwargs)
		self.queue = queue
	@pyqtSlot()
	def run(self):
		while True:
			text = self.queue.get()
			self.mySignal.emit(text)
if __name__ == "__main__":
	queue = Queue()
	sys.stdout = WriteStream(queue)

	qapp = QtGui.QApplication(sys.argv)
	form = MainClass()
	form.show()
	thread = QtCore.QThread()
	reciever = MyReciever(queue)
	reciever.mySignal.connect(form.append_text)
	reciever.moveToThread(thread)
	thread.started.connect(reciever.run)
	thread.start()
	
	qapp.exec_()
