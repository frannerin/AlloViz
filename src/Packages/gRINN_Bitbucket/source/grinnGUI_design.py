# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'grinnGUI.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_gRINN(object):
    def setupUi(self, gRINN):
        gRINN.setObjectName("gRINN")
        gRINN.resize(659, 447)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(gRINN.sizePolicy().hasHeightForWidth())
        gRINN.setSizePolicy(sizePolicy)
        gRINN.setMaximumSize(QtCore.QSize(659, 447))
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("clover.ico"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        gRINN.setWindowIcon(icon)
        self.centralWidget = QtWidgets.QWidget(gRINN)
        self.centralWidget.setObjectName("centralWidget")
        self.label = QtWidgets.QLabel(self.centralWidget)
        self.label.setGeometry(QtCore.QRect(50, 10, 161, 151))
        self.label.setText("")
        self.label.setPixmap(QtGui.QPixmap("clover.ico"))
        self.label.setScaledContents(True)
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(self.centralWidget)
        self.label_2.setGeometry(QtCore.QRect(240, 10, 411, 161))
        font = QtGui.QFont()
        font.setPointSize(36)
        self.label_2.setFont(font)
        self.label_2.setAcceptDrops(True)
        self.label_2.setTextFormat(QtCore.Qt.RichText)
        self.label_2.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)
        self.label_2.setWordWrap(True)
        self.label_2.setObjectName("label_2")
        self.pushButton = QtWidgets.QPushButton(self.centralWidget)
        self.pushButton.setGeometry(QtCore.QRect(60, 200, 241, 141))
        font = QtGui.QFont()
        font.setPointSize(24)
        self.pushButton.setFont(font)
        self.pushButton.setObjectName("pushButton")
        self.pushButton_2 = QtWidgets.QPushButton(self.centralWidget)
        self.pushButton_2.setGeometry(QtCore.QRect(360, 200, 241, 141))
        font = QtGui.QFont()
        font.setPointSize(24)
        self.pushButton_2.setFont(font)
        self.pushButton_2.setObjectName("pushButton_2")
        self.label_3 = QtWidgets.QLabel(self.centralWidget)
        self.label_3.setGeometry(QtCore.QRect(70, 340, 111, 41))
        self.label_3.setAlignment(QtCore.Qt.AlignCenter)
        self.label_3.setObjectName("label_3")
        self.label_4 = QtWidgets.QLabel(self.centralWidget)
        self.label_4.setGeometry(QtCore.QRect(340, 340, 111, 41))
        self.label_4.setAlignment(QtCore.Qt.AlignCenter)
        self.label_4.setObjectName("label_4")
        self.label_5 = QtWidgets.QLabel(self.centralWidget)
        self.label_5.setGeometry(QtCore.QRect(480, 340, 111, 41))
        self.label_5.setAlignment(QtCore.Qt.AlignCenter)
        self.label_5.setObjectName("label_5")
        self.label_6 = QtWidgets.QLabel(self.centralWidget)
        self.label_6.setGeometry(QtCore.QRect(210, 340, 101, 41))
        self.label_6.setAlignment(QtCore.Qt.AlignCenter)
        self.label_6.setObjectName("label_6")
        gRINN.setCentralWidget(self.centralWidget)
        self.menuBar = QtWidgets.QMenuBar(gRINN)
        self.menuBar.setGeometry(QtCore.QRect(0, 0, 659, 22))
        self.menuBar.setObjectName("menuBar")
        gRINN.setMenuBar(self.menuBar)
        self.mainToolBar = QtWidgets.QToolBar(gRINN)
        self.mainToolBar.setObjectName("mainToolBar")
        gRINN.addToolBar(QtCore.Qt.TopToolBarArea, self.mainToolBar)
        self.statusBar = QtWidgets.QStatusBar(gRINN)
        self.statusBar.setObjectName("statusBar")
        gRINN.setStatusBar(self.statusBar)

        self.retranslateUi(gRINN)
        QtCore.QMetaObject.connectSlotsByName(gRINN)

    def retranslateUi(self, gRINN):
        _translate = QtCore.QCoreApplication.translate
        gRINN.setWindowTitle(_translate("gRINN", "gRINN"))
        self.label_2.setText(_translate("gRINN", "<html><head/><body><p>gRINN</p><p><span style=\" font-size:24pt;\">get Residue Interaction eNergies and Networks</span></p></body></html>"))
        self.pushButton.setText(_translate("gRINN", "New\n"
" Calculation"))
        self.pushButton_2.setText(_translate("gRINN", "View\n"
" Results"))
        self.label_3.setText(_translate("gRINN", "<html><head/><body><p><a href=\"http://grinn.readthedocs.io/en/latest/tutorial.html\"><span style=\" font-size:20pt; text-decoration: underline; color:#0000ff;\">Tutorial</span></a></p></body></html>"))
        self.label_4.setText(_translate("gRINN", "<html><head/><body><p><a href=\"http://grinn.readthedocs.io/en/latest/credits.html\"><span style=\" font-size:20pt; text-decoration: underline; color:#0000ff;\">Credits</span></a></p></body></html>"))
        self.label_5.setText(_translate("gRINN", "<html><head/><body><p><a href=\"http://grinn.readthedocs.io/en/latest/contact.html\"><span style=\" font-size:20pt; text-decoration: underline; color:#0000ff;\">Contact</span></a></p></body></html>"))
        self.label_6.setText(_translate("gRINN", "<html><head/><body><p><a href=\"http://grinn.readthedocs.io/en/latest/faq.html\"><span style=\" font-size:20pt; text-decoration: underline; color:#0000ff;\">FAQ</span></a></p></body></html>"))

