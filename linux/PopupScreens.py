#!/usr/bin/python


### Description Popup screen

import subprocess
from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import SIGNAL
from PyQt4.Qt import *

import popup_options
import TnCloneDiagnosisPlot
import os
import numpy as np
from PipelineThread import DiagnosisThreadAnalyse


class DescriptionPopupScreen(QtGui.QWidget):
	
	def __init__(self,desc_type=None):
		self.dialects = {
"sample_info" : """The Sample Information contains two columns<br/><br/>
1. Gene name<br/>
2. number of samples<br/><br/>

Both columns are <b>TAB</b> separated<br/><br/>
<b>CAUTION !!!</b><br/><br/>
The Gene name and number of samples to analyse must be same compared to sorting information.<br/><br/>

Example.<br/><br/>
Let's say out gene name is MYGENE and have 12 samples to analyse<br/><br/>
Then we have<br/><br/>
MYGENE	12<br/><br/>

""",
"sort_info" : """The Sorting Information for TnSuite is described below<br/><br/>
There are 5 columns<br/></br><br/>
First  Column  : Gene name ( We recommend to use 6  chracters )<br/>
Second Column  : Number of Samples to analyse<br/>
Third  Column  : Tn5 specific forward index barcode sequence (8-mers)<br/>
Fourth Column  : Tn5 specific reverse index barcode sequence (8-mers)<br/>
Fifth Column   : The RAW NGS read file name (forward direction)<br/>
Sixth  Column  : The RAW NGS read file name (reverse direction)<br/><br/><br/>

<b>CAUTION !!!</b><br/><br/>
I. All columns must be TAB separated<br/>
II. For fourth, fifth column, file list must be comma separated (See example below)<br/><br/>

Lets say we want to analyse a insert gene called MYGENE<br/>
The number of samples are 12<br/>
Tn5 Specific barcode is GTCCATCA<br/>
Forward side NGS files 1_1.fastq,2_1.fastq,3_1.fastq<br/>
Reverse side NGS files 1_2.fastq,2_2.fastq,3_2.fastq<br/><br/>

Then for this sample, we have format below<br/><br/>

MYGENE\t12\tGTCCATCA\t\tGTCCATCA\t1_1.fastq,2_1.fastq,3_1.fastq\t1_2.fastq,2_2.fastq,3_2.fastq<br/>
""",
"bed_info" : """The BED file ( region file ) cusually contains 4 columns
<br/>
<br/>
In Genomics 4 <b>TAB</b> separated coulmns are usually used<br/><br/>
<br/>I.   CHROMOSOME
<br/>II.  START OF REGION
<br/>III. END OF REGION
<br/>IV.  DESCRIPTION FOR THIS REGION (such as intron, exon etc.)
<br/><br/>
But for us we will use first three columns as described above
<br/><br/>
<br/>I.   GENE
<br/>II.  START OF REGION
<br/>III. END OF REGION
<br/><br/>
<b>CAUTION !!!</b>
<br/>The gene name must be same compared to Sorting information ans Sample information
<br/><br/>
If you do not know how to make that file, just following instruction
<br/><br/>
Let's say your TnSuite installation location is <b>TNSUITE</b>
<br/><br/>
Type following command line at terminal window by pressing (ctrl + alt + T)
<br/><br/>
STEP 1 ( Change directory to your destination )
<br/><br/>
cd <b>TNSUITE</b>
<br/>
cd src
<br/>
<br/><br/>
STEP 2 ( Launch accessory program you need )
<br/><br/>
Suppose your reference files are located under <b>INPUT_REFRENCE</b><br/>
Suppose your output directory is <b>OUTPUT_DIR</b><br/>
Suppose your BED file's prefix ( the file name before extension ) is <b>PREFIX</b><br/>
<br/>Then type below<br/><br/>
python bed_generator.py --ipath <b>INPUT_REFERENCE</b> --opath <b>OUTPUT_DIR</b> --prefix <b>PREFIX</b>

<br/><br/>
Then use PREFIX.bed located file at OUTPUT_DIR
""",
"fasta_info" : """
Reference file must be FASTA format
<br/><br/>
For whome do not know FASTA format, there will be shor description below
<br/><br/>
It usually ends with file extension as fa or fasta (i.e. E_coli.fa or human.fasta)
<br/>
File contents are look like below
<br/><br/>
>some_name_to_describe_this_species additional_description
<br/>DNA sequence(or RNA, Protein)
<br/><br/>
Example below for E.Coli
<br/><br/>
>Ecoli TolC knockout strain
<br/>ATGCCGATGACGATGATGATGATAGTAGATGATGAGATG
<br/>CCGATGAGATAGACGATGAGTGGCCCAGTAGACGATAGA
<br/>AACCCGAGTGAGACCC... (and so on)
<br/><br/>
The sequence field can be either multiple lines or a single line
<br/>See details at https://en.wikipedia.org/wiki/FASTA_format<br/>
<br/>
Here, we have some constraints for fasta file<br/>
<br/>
I. At header field which is starts with '>' notation, this must equal to gene name we used for analysis
<br/>
II. The file name must be gene_name.gene.fa which means we have file extension as gene.fa
<br/>
III. Sequence must be in the linear form
<br/>
Example below
<br/>
Suppose we have to analyse <b>MYGENE</b>
<br/>
Then we have a file with name <b>MYGENE</b>.gene.fa<br/>
And the contents are
<br/><br/>
><b>MYGENE</b>
<br/>ATGAGTAGTAGTAGTCAGATAGACCATC.....
<br/><br/>
<b>NOT</b><br/><br/>
><b>MYGENE</b><br/>
ATGAGTAGT
AGTAGTCAG
ATAGACCAT
C.....
<br/><br/> Which is multiple lined sequence format<br/>
"""
}
		QWidget.__init__(self)
		self.edit = QtGui.QTextEdit()
		self.layout = QtGui.QVBoxLayout()
		self.setLayout(self.layout)
		self.edit.setGeometry(QtCore.QRect(0,0,1000,500))
		self.btn = QPushButton("Ok. I understood",self)
		self.btn.setGeometry(QRect(865,560,120,30))
		self.showDialog(desc_type)
		self.layout.addWidget(self.edit)
		self.layout.addWidget(self.btn)
		self.btn.clicked.connect(self.close)
		
	def showDialog(self,desc_type):
		if desc_type=="bed":
			s = self.dialects["bed_info"]
			self.edit.setReadOnly(True)
			self.edit.setHtml(s)
		elif desc_type == "sort":
			s = self.dialects["sort_info"]
			self.edit.setReadOnly(True)
			self.edit.setHtml(s)
		elif desc_type == "sample":
			s = self.dialects["sample_info"]
			self.edit.setReadOnly(True)
			self.edit.setHtml(s)
		elif desc_type == "fasta":
			s = self.dialects["fasta_info"]
			self.edit.setReadOnly(True)
			self.edit.setHtml(s)
### Diagnosis plot

class DiagnosisPlot(QtGui.QWidget, TnCloneDiagnosisPlot.Ui_Form):
	def __init__(self, sample_info_file, trim_path, ref_path):
		super(self.__class__,self).__init__()
		self.sample_info = sample_info_file
		self.trim_path = trim_path
		self.ref_path = ref_path
		self.setupUi(self)
		self.pushButton.clicked.connect(self.diagnosis)
		self.pushButton_2.clicked.connect(self.draw)
		self.pushButton_4.clicked.connect(self.prev)
		self.pushButton_3.clicked.connect(self.next_plot)
		
		self.log_screen = self.textEdit
		self.log_prefix = "[DiagnosisPlot] "
		self.root_path = "/".join(trim_path.split("/")[:-1])
		self.dataset = None
		self.first_draw = False
		self.filenames = None
		self.plot_index = 0
		self.filename_desc = self.lineEdit
	def diagnosis(self):
		self.analyser = DiagnosisThreadAnalyse(self.sample_info ,self.trim_path, self.ref_path,self.log_screen)
		self.connect(self.analyser,SIGNAL('finished()'),self.done)
		self.analyser.start()
	def done(self):
		QtGui.QMessageBox.information(self, "Done!" , "Done Diagnosis!!")
	def _parse_all_pileups(self) :
		### return all parsed data
		pileup_path = "/".join(self.trim_path.split("/")[:-1]) + "/diagnosis/pileup"
		files = os.listdir(pileup_path)
		dataset = {}
		for filename in files :
			self.log_screen.append(self.log_prefix + " Pileup-file : %s"%(filename))
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
	def _sub_plot_command(self,dataset,filenames):
		self.log_screen.append(self.log_prefix + "Preparing data for drawing figure...")
		canvas_ax = self.widget.canvas.ax
		canvas = self.widget.canvas
		target_filename = filenames[self.plot_index]
		self.filename_desc.setText(target_filename)
		canvas_ax.clear()
		posits = dataset[ filenames[self.plot_index] ]["posits"]
		depths = dataset[ filenames[self.plot_index] ]["depths"]
		xmin = min(posits)
		xmax = max(posits)
		canvas_ax.set_ylabel('Coverage')
		canvas_ax.set_xlabel('Coordinates')
		#import matplotlib.pyplot as plt
		#axes = plt.gca()
		canvas_ax.set_xlim([xmin,xmax])
		#axes.set_xlim([xmin, xmax])
		canvas_ax.plot(posits, depths,lw=2)
		canvas_ax.axhline(y=np.mean(depths), ls='dashed', color='red', lw = 3)
		canvas_ax.fill_between(posits, depths, facecolor='blue', alpha=0.5)
		#for label in canvas_ax.get_yticklabels():
		#	label.set_visible(False)
		self.log_screen.append(self.log_prefix + "Done preparing data for drawing figure...")
		canvas_ax.set_title('Diagnostic plot for evaluating coverage of given reference\n(Dashed line indicates average coverage)')
		canvas.draw()

	def draw(self):
		
		if self.first_draw == False:
			self.log_screen.append(self.log_prefix + "Parsing all data from previously diagnosis")
			self._parse_all_pileups()
			dataset = self.dataset

			filenames = dataset.keys()
			filenames.sort()
			self.filenames = filenames
			### see _sub_plot_command function
			self._sub_plot_command(dataset,filenames)
			self.first_draw = True
		else:
			dataset = self.dataset

			filenames = self.filenames

			self._sub_plot_command(dataset,filenames)

	def prev(self):
		self.log_screen.append(self.log_prefix + "Move to previous figure...")
		if self.first_draw == False:
			self.draw()
			self.plot_index -= 1
		else:
			filenames = self.filenames
			maxrange = len(filenames)

			if self.plot_index == 0 :
				self.plot_index = maxrange - 1
				self.draw()
				self.plot_index -= 1
			else:
				self.draw()
				self.plot_index -= 1

	def next_plot(self):
		self.log_screen.append(self.log_prefix + "Move to next figure...")
		if self.first_draw == False:
			self.draw()
			self.plot_index += 1
		else:
			if self.filenames != None:
				filenames = self.filenames
			else:
				self.filenames = dataset.keys()
				self.filenames.sort()
				filenames = self.filenames
			maxrange = len(filenames)
			if self.plot_index == maxrange-1 :
				self.plot_index = 0
				self.draw()
				self.plot_index += 1
			else:
				self.draw()
				self.plot_index += 1

### Option popup screen
class OptionPopupScreen(QtGui.QWidget, popup_options.Ui_Form):
	""" Required field here """
	### ASSEMBLY
	kmer_update = QtCore.pyqtSignal(str)
	mink_occ = QtCore.pyqtSignal(str)
	sort_on = QtCore.pyqtSignal(bool)
	trim_on = QtCore.pyqtSignal(bool)
	analysis_on = QtCore.pyqtSignal(bool)
	assembly_on = QtCore.pyqtSignal(bool)
	bed_path = QtCore.pyqtSignal(str)
	### ALIGNMENT
	match_score = QtCore.pyqtSignal(str) 
	mismatch_score = QtCore.pyqtSignal(str)
	gap_open = QtCore.pyqtSignal(str)
	gap_extension = QtCore.pyqtSignal(str)

	### Other misc options
	down_sample = QtCore.pyqtSignal(bool)
	down_sample_rate = QtCore.pyqtSignal(float)
	depth_correction = QtCore.pyqtSignal(bool)
	mismatch_correction = QtCore.pyqtSignal(bool)
	max_mis = QtCore.pyqtSignal(int)
	gpu_avail = QtCore.pyqtSignal(bool)
	def __init__(self):
		super(self.__class__,self).__init__()
		self.setupUi(self)
		### Browse button
		self.pushButton_4.clicked.connect(self.browse_bedfile)
		### Save changes
		self.pushButton_10.clicked.connect(self.save_changes)
		self.pushButton_11.clicked.connect(self.set_default)
		self.pushButton.clicked.connect(lambda:self.description('kmer'))
		self.pushButton_2.clicked.connect(lambda:self.description('mink'))
		self.pushButton_3.clicked.connect(lambda:self.description('step'))
		self.pushButton_5.clicked.connect(lambda:self.description('bed'))
		self.pushButton_6.clicked.connect(lambda:self.description('match'))
		self.pushButton_7.clicked.connect(lambda:self.description('mismatch'))
		self.pushButton_8.clicked.connect(lambda:self.description('go'))
		self.pushButton_9.clicked.connect(lambda:self.description('ge'))
		self.screen = None
		self.bedScreen = None
	def browse_bedfile(self):
		bedpath = QtGui.QFileDialog.getOpenFileName(self,"Browse BED file", "/" , "*")
		self.lineEdit_17.setText(bedpath)
	def description(self, dtype):
		self.screen = DescriptionPopupScreen()
		if dtype == 'kmer' :
			s = "The length of k-mer. If you don't use bed file, then start/end sequence must be equal to this value"
			self.screen.edit.setReadOnly(True)
			self.screen.edit.setText(s)
		elif dtype == 'mink' :
			s = "The minimum k-mer occurence. If the count of k-mer are smaller than this value, those will be discarded."
			self.screen.edit.setReadOnly(True)
			self.screen.edit.setText(s)
		elif dtype =='step':
			s = "The analysis step. Check what you want. But proccess are streamlines so don't uncheck box in the middle"
			self.screen.edit.setReadOnly(True)
			self.screen.edit.setText(s)
		elif dtype == 'bed':
			s = self.screen.dialects['bed_info']
			
			
			self.screen.setGeometry(QRect(0,0,1000,600))
			self.screen.edit.setReadOnly(True)
			self.screen.edit.setHtml(s)

		elif dtype == 'match' :
			s = 'Alignement match score.\nValue is integer.'
			self.screen.edit.setReadOnly(True)
			self.screen.edit.setText(s)
		elif dtype == 'mismatch' :
			s = 'Alignement mis-match score.\nValue is integer.'
			self.screen.edit.setReadOnly(True)
			self.screen.edit.setText(s)
		elif dtype == 'go' :
			s = 'Alignement gap open penalty.\nValue is integer.'
			self.screen.edit.setReadOnly(True)
			self.screen.edit.setText(s)
		elif dtype == 'ge' :
			s = 'Alignement gap extension penalty.\nValue is integer.'
			self.screen.edit.setReadOnly(True)
			self.screen.edit.setText(s)
		self.screen.show()
	def set_default(self):
		self.lineEdit_3.setText('63')
		self.lineEdit_5.setText('3')
		self.checkBox.setChecked(True)
		self.checkBox_2.setChecked(True)
		self.checkBox_3.setChecked(True)
		self.checkBox_4.setChecked(True)
		self.lineEdit_17.setText('')
		self.lineEdit_9.setText('1')
		self.lineEdit_11.setText('3')
		self.lineEdit_13.setText('5')
		self.lineEdit_15.setText('3')
	def save_changes(self):
		choice = QtGui.QMessageBox.question(self,"Save Changes", "Are you sure to save changes and exit?", QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
		if choice == QtGui.QMessageBox.Yes:
			self.change_vals()
			self.close()
		else:
			self.clear_updates()
			self.close()
	def change_vals(self):
		kmer_update_val = str(self.lineEdit_3.text())
		mink_occ_val = str(self.lineEdit_5.text())
		sort_on_val = bool(self.checkBox.isChecked())
		trim_on_val = bool(self.checkBox_2.isChecked())
		assembly_on_val = bool(self.checkBox_3.isChecked())
		analysis_on_val = bool(self.checkBox_4.isChecked())
		bed_path_val = str(self.lineEdit_17.text())
		match_score_val = str(self.lineEdit_9.text())
		mismatch_score_val = str(self.lineEdit_11.text())
		gap_open_val = str(self.lineEdit_13.text())
		gap_extension_val = str(self.lineEdit_15.text())

		down_sample = bool(self.checkBox_5.isChecked())
		down_sample_rate = float(self.lineEdit_20.text())
		depth_correction = bool(self.checkBox_6.isChecked())
		mismatch_correction = bool(self.checkBox_7.isChecked())
		max_mis = int(self.lineEdit_23.text())
		#gpu_avail = bool(self.checkBox_8.isChecked())

		self.kmer_update.emit(kmer_update_val)
		self.mink_occ.emit(mink_occ_val)
		self.sort_on.emit(sort_on_val)
		self.trim_on.emit(trim_on_val)
		self.assembly_on.emit(assembly_on_val)
		self.analysis_on.emit(analysis_on_val)
		self.bed_path.emit(bed_path_val)
		self.match_score.emit(match_score_val)
		self.mismatch_score.emit(mismatch_score_val)
		self.gap_open.emit(gap_open_val)
		self.gap_extension.emit(gap_extension_val)
		self.down_sample.emit(down_sample)
		self.down_sample_rate.emit(down_sample_rate)
		self.depth_correction.emit(depth_correction)
		self.mismatch_correction.emit(mismatch_correction)
		self.max_mis.emit(max_mis)
		#self.gpu_avail.emit(gpu_avail)

	def clear_updates(self):
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
