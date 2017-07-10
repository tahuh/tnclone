# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'TnCloneDiagnosisPlot.ui'
#
# Created: Wed Nov 23 17:13:17 2016
#      by: PyQt4 UI code generator 4.10.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName(_fromUtf8("Form"))
        Form.resize(640, 630)
        self.pushButton = QtGui.QPushButton(Form)
        self.pushButton.setGeometry(QtCore.QRect(10, 10, 141, 27))
        self.pushButton.setObjectName(_fromUtf8("pushButton"))
        self.pushButton_2 = QtGui.QPushButton(Form)
        self.pushButton_2.setGeometry(QtCore.QRect(160, 10, 211, 27))
        self.pushButton_2.setObjectName(_fromUtf8("pushButton_2"))
        self.pushButton_3 = QtGui.QPushButton(Form)
        self.pushButton_3.setGeometry(QtCore.QRect(510, 10, 98, 27))
        self.pushButton_3.setObjectName(_fromUtf8("pushButton_3"))
        self.pushButton_4 = QtGui.QPushButton(Form)
        self.pushButton_4.setGeometry(QtCore.QRect(410, 10, 98, 27))
        self.pushButton_4.setObjectName(_fromUtf8("pushButton_4"))
        self.widget = matplotlibWidget(Form)
        self.widget.setGeometry(QtCore.QRect(10, 70, 601, 421))
        self.widget.setObjectName(_fromUtf8("widget"))
        self.textEdit = QtGui.QTextEdit(Form)
        self.textEdit.setGeometry(QtCore.QRect(10, 500, 601, 121))
        self.textEdit.setObjectName(_fromUtf8("textEdit"))
        self.lineEdit = QtGui.QLineEdit(Form)
        self.lineEdit.setGeometry(QtCore.QRect(160, 40, 451, 27))
        self.lineEdit.setObjectName(_fromUtf8("lineEdit"))
        self.lineEdit_2 = QtGui.QLineEdit(Form)
        self.lineEdit_2.setGeometry(QtCore.QRect(10, 40, 141, 27))
        self.lineEdit_2.setObjectName(_fromUtf8("lineEdit_2"))

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "TnClone-Diagnosis Plot", None))
        self.pushButton.setText(_translate("Form", "Diagnosis samples", None))
        self.pushButton_2.setText(_translate("Form", "Plot First Sample", None))
        self.pushButton_3.setText(_translate("Form", "Next", None))
        self.pushButton_4.setText(_translate("Form", "Previous", None))
        self.lineEdit_2.setText(_translate("Form", "FILE processed", None))

from matplotlibwidget import matplotlibWidget
