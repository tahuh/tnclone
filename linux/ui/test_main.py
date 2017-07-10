#!/usr/bin/python

import sys
from PyQt4 import QtGui , QtCore
from PyQt4.QtCore import SIGNAL
from PyQt4.Qt import *
import design_test


class PipelineThread(QtCore.QThread):
	def __init__(self, m, M) :
		QtCore.QThread.__init__(self)
		self.m = m
		self.M = M

	def __del__(self):
		self.wait()

	def run(self) :
		s = 0
		for x in xrange(self.m, self.M):
			s += x

### Pupup-window

class MyPopUpSortingInformation(QWidget):
	def __init__(self):
		QWidget.__init__(self)
		self.edit = QtGui.QTextEdit()
		self.layout = QtGui.QVBoxLayout()
		self.setLayout(self.layout)
		self.edit.setGeometry(QtCore.QRect(0,0,1000,500))
		self.btn = QPushButton("Ok. I understood",self)
		self.btn.setGeometry(QRect(865,560,120,30))
		#self.btn.setStyleSheet("background-color: brown")
		self.showDialog()
		self.layout.addWidget(self.edit)
		self.layout.addWidget(self.btn)
		self.btn.clicked.connect(self.close)

	def showDialog(self):

		sentences = """The Sorting Information for TnSuite is described below<br/><br/>
There are 5 columns<br/></br><br/>
First  Column  : Gene name ( We recommend to use 6  chracters )<br/>
Second Column  : Number of Samples to analyse<br/>
Third  Column  : Tn5 specific barcode (8-mers)<br/>
Fourth Column  : The RAW NGS read file name (forward direction)<br/>
Fifth  Column  : The RAW NGS read file name (reverse direction)<br/><br/><br/>

<b>CAUTION !!!</b><br/><br/>
I. All columns must be TAB separated<br/>
II. For fourth, fifth column, file list must be comma separated (See example below)<br/><br/>

Lets say we want to analyse a insert gene called MYGENE<br/>
The number of samples are 12<br/>
Tn5 Specific barcode is GTCCATCA<br/>
Forward side NGS files 1_1.fastq,2_1.fastq,3_1.fastq<br/>
Reverse side NGS files 1_2.fastq,2_2.fastq,3_2.fastq<br/><br/>

Then for this sample, we have format below<br/><br/>

MYGENE\t12\tGTCCATCA\t1_1.fastq,2_1.fastq,3_1.fastq\t1_2.fastq,2_2.fastq,3_2.fastq<br/>
""".split("\n")
		S = "\n".join(sentences)
		self.edit.setReadOnly(True)
		#self.edit.setPlainText(S)
		self.edit.setHtml(S)
		#for s in sentences :
		#	print s
		#	self.edit.append(s)
		
class MyPopUpSampleInformation(QWidget):
	def __init__(self):
		QWidget.__init__(self)
		self.edit = QtGui.QTextEdit()
		self.layout = QtGui.QVBoxLayout()
		self.setLayout(self.layout)
		self.edit.setGeometry(QtCore.QRect(0,0,1000,500))
		self.btn = QPushButton("Ok. I understood",self)
		self.btn.setGeometry(QRect(865,560,120,30))
		self.showDialog()
		self.layout.addWidget(self.edit)
		self.layout.addWidget(self.btn)
		self.btn.clicked.connect(self.close)

	def showDialog(self):

		sentences = """The Sample Information contains two columns<br/><br/>
1. Gene name<br/>
2. number of samples<br/><br/>

Both columns are <b>TAB</b> separated<br/><br/>
<b>CAUTION !!!</b><br/><br/>
The Gene name and number of samples to analyse must be same compared to sorting information.<br/><br/>

Example.<br/><br/>
Let's say out gene name is MYGENE and have 12 samples to analyse<br/><br/>
Then we have<br/><br/>
MYGENE	12<br/><br/>

""".split("\n")
		S = "\n".join(sentences)
		self.edit.setReadOnly(True)
#		self.edit.setPlainText(S)
		self.edit.setHtml(S)

class MyPopUpBedInformation(QWidget):
	def __init__(self):
		QWidget.__init__(self)
		self.edit = QtGui.QTextEdit()
		self.layout = QtGui.QVBoxLayout()
		self.setLayout(self.layout)
		self.edit.setGeometry(QtCore.QRect(0,0,1000,500))
		self.btn = QPushButton("Ok. I understood",self)
		self.btn.setGeometry(QRect(865,560,120,30))
		self.showDialog()
		self.layout.addWidget(self.edit)
		self.layout.addWidget(self.btn)
		self.btn.clicked.connect(self.close)

	def showDialog(self):

		sentences = """
The BED file ( region file ) cusually contains 4 columns
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
""".split("\n")
		S = "\n".join(sentences)
		self.edit.setReadOnly(True)
#		self.edit.setPlainText(S)
		self.edit.setHtml(S)

class MyPopUpRefInformation(QWidget):
	def __init__(self):
		QWidget.__init__(self)
		self.edit = QtGui.QTextEdit()
		self.layout = QtGui.QVBoxLayout()
		self.setLayout(self.layout)
		self.edit.setGeometry(QtCore.QRect(0,0,1000,500))
		self.btn = QPushButton("Ok. I understood",self)
		self.btn.setGeometry(QRect(865,560,120,30))
		self.showDialog()
		self.layout.addWidget(self.edit)
		self.layout.addWidget(self.btn)
		self.btn.clicked.connect(self.close)

	def showDialog(self):

		sentences = """
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
""".split("\n")
		S = "\n".join(sentences)
		self.edit.setReadOnly(True)
#		self.edit.setPlainText(S)
		self.edit.setHtml(S)

import untitled
### A Class shows options we need

class ShowOptionClass(QtGui.QWidget, untitled.Ui_Form):
	def __init__(self):
		super(self.__class__,self).__init__()
		self.setupUi(self)
#		QWidget.__init__(self)
		
		#self.edit = QtGui.QTextEdit()
#		self.layout = QtGui.QVBoxLayout()
#		self.layout_main = QtGui.QHBoxLayout()
#		self.setLayout(self.layout_main)
#		self.layout_main.setLayout(self.layout)
		#self.edit.setGeometry(QtCore.QRect(0,0,1000,500))
#		self.btn = QPushButton("Add these option field",self)
#		self.btn.setGeometry(QRect(865,560,120,30))
#		self.btn2 = QPushButton("What the fuck!")
#		self.btn2.setGeometry(QRect(765,460,120,30))
		#self.showDialog()
#		#self.layout.addWidget(self.edit)
#		self.layout.addWidget(self.btn2)
#		self.layout.addWidget(self.btn)

#		self.btn.clicked.connect(self.close)

	
class Body(QtGui.QMainWindow, design_test.Ui_MainWindow):
	def __init__(self):
		super(self.__class__,self).__init__()
		self.setupUi(self)

		self.pushButton_15.clicked.connect(self.abortapp)
		self.pushButton_7.clicked.connect(self.launch)
		self.pushButton_8.clicked.connect(self.quitapp)

		self.pushButton_12.clicked.connect(self.popSort)
		self.pushButton_13.clicked.connect(self.popSample)
		self.pushButton_14.clicked.connect(self.popBed)
		self.pushButton_16.clicked.connect(self.popRef)
		self.pushButton_10.clicked.connect(self.popOptionDialog)
		self.w = None
		self.w2 = None
		self.w3 = None
		self.w4 = None
		self.optionDialog = None
	def popOptionDialog(self):
		self.optionDialog = ShowOptionClass()
		self.optionDialog.setWindowTitle("Additional Options....")
		#self.optionDialog.setGeometry(QRect(0,0,1000,600))
		self.optionDialog.show()
	def popSort(self):
		self.w = MyPopUpSortingInformation()
		self.w.setWindowTitle("Sorting information file format description")
		self.w.setGeometry(QRect(0,0,1000,600))
		self.w.showDialog()
		self.w.show()
	def popSample(self):
		self.w2 = MyPopUpSampleInformation()
		self.w2.setWindowTitle("Sample information file format description")
		self.w2.setGeometry(QRect(0,0,1000,600))
		self.w2.showDialog()
		self.w2.show()

	def popBed(self):
		self.w3 = MyPopUpBedInformation()
		self.w3.setWindowTitle("Bed information file format description")
		self.w3.setGeometry(QRect(0,0,1000,600))
		self.w3.showDialog()
		self.w3.show()

	def popRef(self):
		self.w4 = MyPopUpRefInformation()
		self.w4.setWindowTitle("Reference information file format description")
		self.w4.setGeometry(QRect(0,0,1000,600))
		self.w4.showDialog()
		self.w4.show()
	def launch(self):
		self.pipe_thread = PipelineThread(10,10000)
		self.connect(self.pipe_thread, SIGNAL("finished()"), self.done)
		self.pipe_thread.start()
		
	def done(self) :
		print "Done!"

	def quitapp(self) :
		choice = QtGui.QMessageBox.question(self.centralwidget, "Exit Prgram", "Are you sure to exit the program?" , QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
		if choice == QtGui.QMessageBox.Yes:
			sys.exit()
		else:
			pass

	def abortapp(self):
		choice = QtGui.QMessageBox.question(self.centralwidget, "Abort program", "Are you sure to abort?" , QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)

		if choice == QtGui.QMessageBox.Yes:
#			self.pipe_thread.stop()
			self.pipe_thread.quit()
			self.pipe_thread.wait()
			print "Aborting..."
			QtGui.QMessageBox.information(self, "Abort!", "Aborted...")
		else:
			pass
if __name__ == "__main__":
	app = QtGui.QApplication(sys.argv)
	form = Body()
	form.show()
	app.exec_()	
