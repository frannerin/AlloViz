#!/usr/bin/env python
# -*- coding: utf-8 -*-
from PyQt5 import QtGui
from PyQt5 import QtCore
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QMessageBox
#import resultsGUI, calcGUI, grinnGUI_design, # Changed by frannerin
from . import calc, corr, common, grinnGUI_design # Changed by frannerin
import sys, time, os, argparse, multiprocessing, subprocess
multiprocessing.freeze_support() # Required for Windows compatibility, harmless for Unix.

class DesignInteract(QMainWindow,grinnGUI_design.Ui_gRINN):

	def __init__(self,parent=None):
		super(DesignInteract,self).__init__(parent)
		self.setupUi(self)

		_translate = QtCore.QCoreApplication.translate
		self.label_3.setText(_translate("gRINN", "<html><head/><body><p><a href=\"http://grinn.readthedocs.io/en/latest/tutorial.html\"><span style=\" font-size:20pt; text-decoration: underline; color:#0000ff;\">Tutorial</span></a></p></body></html>"))
		self.label_4.setText(_translate("gRINN", "<html><head/><body><p><a href=\"http://grinn.readthedocs.io/en/latest/credits.html\"><span style=\" font-size:20pt; text-decoration: underline; color:#0000ff;\">Credits</span></a></p></body></html>"))
		self.label_5.setText(_translate("gRINN", "<html><head/><body><p><a href=\"http://grinn.readthedocs.io/en/latest/contact.html\"><span style=\" font-size:20pt; text-decoration: underline; color:#0000ff;\">Contact</span></a></p></body></html>"))
		self.label_6.setText(_translate("gRINN", "<html><head/><body><p><a href=\"http://grinn.readthedocs.io/en/latest/faq.html\"><span style=\" font-size:20pt; text-decoration: underline; color:#0000ff;\">FAQ</span></a></p></body></html>"))

		self.label_3.setOpenExternalLinks(True)
		self.label_4.setOpenExternalLinks(True)
		self.label_5.setOpenExternalLinks(True)
		self.label_6.setOpenExternalLinks(True)
		self.pushButton.clicked.connect(self.calculateGUI)
		self.pushButton_2.clicked.connect(self.resultsGUI)

	def calculateGUI(self):
		self.formGetResIntEnGUI = calcGUI.DesignInteractCalculate(self)
		self.formGetResIntEnGUI.show()
		icon = QtGui.QIcon()
		pixmap = QtGui.QPixmap(common.resource_path(
			os.path.join('resources','clover.ico')))
		icon.addPixmap(pixmap,QtGui.QIcon.Normal, QtGui.QIcon.Off)
		self.formGetResIntEnGUI.setWindowIcon(icon)
		self.formGetResIntEnGUI.label_3.setPixmap(pixmap)

	def resultsGUI(self):
		self.formResults = resultsGUI.DesignInteractResults(self)
		self.formResults.show()
		icon = QtGui.QIcon()
		pixmap = QtGui.QPixmap(common.resource_path(
			os.path.join('resources','clover.ico')))
		icon.addPixmap(pixmap,QtGui.QIcon.Normal, QtGui.QIcon.Off)
		self.formResults.setWindowIcon(icon)
		# Skip through tab widgets to show each GUI component 
		# (apparently necessary for plots to draw correctly...
		self.formResults.tabWidget.setCurrentIndex(0)
		self.formResults.tabWidget.setCurrentIndex(2)
		self.formResults.tabWidget.setCurrentIndex(3)
		self.formResults.tabWidget.setCurrentIndex(4)
		self.formResults.tabWidget.setCurrentIndex(5)
		self.formResults.tabWidget_2.setCurrentIndex(0)
		self.formResults.tabWidget_2.setCurrentIndex(1)
		self.formResults.tabWidget_2.setCurrentIndex(0)
		self.formResults.tabWidget.setCurrentIndex(0)
		time.sleep(1)
		folderLoaded = self.formResults.updateOutputFolder()
		#if not folderLoaded:
		#	self.formResults.close()

	def closeEvent(self, event):
		message = False
		closeCalc = False
		# Check whether any calcGUI or resultsGUI views have been created.
		if hasattr(self,"formGetResIntEnGUI"):
			if self.formGetResIntEnGUI.isVisible():
				message = 'At least one "New Calculation" interface is active. Are you sure ?'
				closeCalc = True
		if hasattr(self,"formResults"):
			if self.formResults.isVisible():
				message = 'At least one "View Results" interface is active. Are you sure ?'

		if message:
			# Is the user sure about this?
			buttonReply = QMessageBox.question(
				self, 'Are you sure?', message,
				QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
			if buttonReply == QMessageBox.No:
				event.ignore()
			elif buttonReply == QMessageBox.Yes:
				if closeCalc:
					self.formGetResIntEnGUI.close()
					event.accept()
				else:
					event.accept()
		elif not message:
			event.accept()

def prepareEnvironment():
	# Set some environment variable for pyinstaller executable function.
	filePath = os.path.abspath(__file__)
	os.environ['FONTCONFIG_FILE'] = common.resource_path(
		os.path.join('data','etc','fonts','fonts.conf'))
	os.environ['FONTCONFIG_PATH'] = common.resource_path(
		os.path.join('data','etc','fonts'))
	os.environ['QT_XKB_CONFIG_ROOT'] = common.resource_path(
		os.path.join('data','xkb'))
	if os.getenv('LD_LIBRARY_PATH'):
		os.environ['LD_LIBRARY_PATH'] = os.getenv('LD_LIBRARY_PATH')+':'+common.resource_path(
		os.path.join('data','xcbglintegrations'))
	else:
		os.environ['LD_LIBRARY_PATH'] = common.resource_path(
		os.path.join('data','xcbglintegrations'))

	os.environ['LD_LIBRARY_PATH'] = os.getenv('LD_LIBRARY_PATH')+':'+common.resource_path(
		'data')
	os.environ['LD_LIBRARY_PATH'] = os.getenv('LD_LIBRARY_PATH')+':'+common.resource_path(
		os.path.join('platforms'))
	if os.getenv('QT_QPA_PLATFORM_PLUGIN_PATH'):
		os.environ['QT_QPA_PLATFORM_PLUGIN_PATH'] = \
		os.getenv('QT_QPA_PLATFORM_PLUGIN_PATH')+':'+common.resource_path('plugins')
	else:
		os.environ['QT_QPA_PLATFORM_PLUGIN_PATH'] = common.resource_path('plugins')
	print('LD_LIBRARY_PATH is set to:')
	print(os.getenv('LD_LIBRARY_PATH'))
	print('QT_QPA_PLATFORM_PLUGIN_PATH is set to:')
	print(os.getenv('QT_QPA_PLATFORM_PLUGIN_PATH'))

def main():
	sys_argv = sys.argv
	sys_argv += ['--style','Fusion']
	app = QApplication(sys_argv)
	form = DesignInteract()
	icon = QtGui.QIcon()
	pixmap = QtGui.QPixmap(common.resource_path(
		os.path.join('resources','clover.ico')))
	icon.addPixmap(pixmap,QtGui.QIcon.Normal, QtGui.QIcon.Off)
	form.setWindowIcon(icon)
	form.label.setGeometry(QtCore.QRect(50, 10, 161, 151))
	form.label.setText("")
	form.label.setScaledContents(True)
	form.label.setPixmap(pixmap)
	form.show()
	app.exec_()

def convert_arg_line_to_args(arg_line):
	# To override the same method of the ArgumentParser (to read options from a file)
	for arg in arg_line.split():
		if not arg.strip():
			continue
		yield arg

# if __name__ == '__main__': # Changed by frannerin
def arg_parser(argsl=None): # Changed by frannerin

	# Construct an argument parser.
	parser = argparse.ArgumentParser(fromfile_prefix_chars='@',
		description='gRINN: get Residue Interaction eNergies and Networks. '
		'gRINN calculates pairwise molecular-mechanics interaction energies '
		'between amino acid residues in the context of molecular dynamics '
		'simulations. gRINN also calculates equal-time correlations between '
		'interaction energies and constructs protein energy networks. '
		'gRINN offers separate graphical user interfaces for calculation of '
		'energies and visualization of results.')

	# Overriding convert_arg_line_to_args in the input parser.
	parser.convert_arg_line_to_args = convert_arg_line_to_args

	parser.add_argument('-calc',action='store_true',default=False,
		help='Activates interaction energy calculation mode. ')

	parser.add_argument('-corr',action='store_true',default=False,
		help='Activates interaction energy correlation calculation mode. ')

	parser.add_argument('-results',action='store_true',default=False,
		help='Activates viewing of results mode. ')

	parser.add_argument('--pdb',default=[False],type=str,nargs=1,
		help='Name of the corresponding PDB file of the DCD trajectory. '
		'Applies only to NAMD-generated trajectories.')

	parser.add_argument('--tpr',default=[False],type=str,nargs=1,
		help='Name of the corresponding TPR file of the XTC/TRR trajectory. '
		'Applies only to GROMACS-generated trajectories.')

	parser.add_argument('--top',default=[False],type=str,nargs=1,
		help='Name of the corresponding PSF file of the DCD trajectory or '
		'TOP file of the XTC/TRR trajectory')

	parser.add_argument('--traj',default=[False],type=str,nargs=1,
		help='Name of the trajectory file')

	parser.add_argument('--numcores',default=[multiprocessing.cpu_count()],
		type=int,nargs=1,
		help='Number of CPU cores to be employed for energy calculation. '
		'If not specified, it defaults to the number of cpu cores present '
		'in your computer.')

	parser.add_argument('--dielectric',default=[1],type=int,nargs=1,
		help='Solute dielectric constant to be used in electrostatic interaction energy '
		'computation. Applies only to NAMD-generated trajectories (DCD).')

	parser.add_argument('--sel1',default=['all'],nargs='+',
		help='A ProDy atom selection string which determines the first group of selected '
		'residues. ')

	parser.add_argument('--sel2',default=['all'],type=str,nargs='+',
		help='A ProDy atom selection string which determines the second group of selected '
		'residues.')

	parser.add_argument('--pairfiltercutoff',type=float,default=[15],nargs=1,
		help='Cutoff distance (angstroms) for pairwise interaction energy calculations. '
		'If not specified, it defaults to 15 Angstroms. '
		'Only those residues that are within the PAIRFILTERCUTOFF distance of each other '
		'for at least PAIRCUTOFFPERCENTAGE of the trajectory will be included '
		'in energy calculations.')

	parser.add_argument('--pairfilterpercentage',default=[75],nargs=1,
		help='When given, residues that	are within the PAIRFILTERCUTOFF distance from each '
		'other for at least PAIRFILTERPERCENTAGE percent of the trajectory will be taken '
		'into account in further evaluations. When not given, it defaults to 75%%)')

	parser.add_argument('--cutoff',default=[12],nargs=1,
		help='Non-bonded interaction distance cutoff (Angstroms) for NAMD-type data.')

	parser.add_argument('--stride',default=[1],type=int,nargs=1,
		help='If specified, a stride with value of STRIDE will be applied to the trajectory '
		'during interaction energy calculation.')

	parser.add_argument('--framerange',type=int,default=[False],nargs='+',
		help='If specified, then only FRAMERANGE\ section of the trajectory will be handled. '
		'For example, if you specify --framerange 100 1000, then only frames between 100 '
		' and 1000 will be included in all calculations. Applies only to grinn -calc '
		'<arguments> calls.')

	parser.add_argument('--exe',default=[None],type=str,nargs=1,
		help='Path to the namd2/gmx executable. Defaults to namd2 or gmx, depending on '
		'whether you specify NAMD2 or GMX type trajectory data (assumes namd2 or gmx is '
		'in the executable search path.')

	parser.add_argument('--parameterfile',default=[False],type=str,nargs='+',
		help='Path to the parameter file(s). Applies only to NAMD-generated data. '
		'You can specify multiple parameters one after each other by placing a blank '
		'space between them, e.g. --parameterfile file1.inp file2.inp. Applies only to '
		'grinn -calc <arguments> calls.')

	parser.add_argument('--calccorr',action='store_true',default=False,
		help='When specified, interaction energy correlation is also calculated following '
		'interaction energy calculation in -calc mode. Equivalent to a grinn -corr call '
		'after grinn -calc. Applies only to grinn -calc <arguments> calls.')

	parser.add_argument('--corrinfile',type=str,nargs=1,help='Path to the CSV file where interaction\
		energies are located in')

	parser.add_argument('--corrintencutoff',default=[1],type=float,nargs=1,help='\
		Mean (average) interaction energy cutoff for filtering interaction energies '
		'(kcal/mol) prior to correlation calculation. If an interaction energy time series '
		'absolute average value is below this cutoff, that interaction energy will not be '
		'taken in correlation calculations.	Defaults to 1 kcal/mol. Applied to grinn -calc '
		'<arguments> and grinn -corr <arguments> calls.')

	parser.add_argument('--corrprefix',type=str,nargs=1,default=[''],
		help='Prefix to the file names for storing calculation results.')

	parser.add_argument('--outfolder',default=[os.path.join(os.getcwd(),
		'grinn_output')],type=str,nargs=1,
		help='Folder path for storing calculation results. If not specified, a folder named '
		'grinn_output will be created in the current working folder. Applies only to grinn '
		'-calc <arguments> calls.')

	parser.add_argument('--version',action='store_true',default=False,
		help='Prints the version number.')

	# Parse arguments.
	args = parser.parse_args(argsl) # Changed by frannerin

	return args # Changed by frannerin

if __name__ == '__main__': # Changed by frannerin
	args = arg_parser(sys.argv) # Changed by frannerin

	calcMode = args.calc
	corrMode = args.corr
	resultsMode = args.results

	# Check which mode is requested. Either one is selected or none is selected to
	# enter GUI mode.
	if [calcMode,corrMode,resultsMode].count(True) > 1:
		print('You should specify either -calc or -corr or specify none of them '
			'to enter the GUI mode.')
		sys.exit(0)
	elif [calcMode,corrMode,resultsMode].count(True) == 1:
		if calcMode:
			# User requested command-line calculation of interaction energies.
			calc.getResIntEn(args)
		elif corrMode:
			# User requested command-line calculation of interaction energy correlations.
			corr.getResIntCorr(args,logFile=None)
		elif resultsMode:
			# User requested to view results of a completed calculation.
			prepareEnvironment()
			resultsGUI.main()
	else:
		if args.version:
			# User requested printing of version.
			versionfile = open(common.resource_path('VERSION'),'r')
			versionline = versionfile.readlines()
			version = versionline[0].rstrip('\n')
			print(version)
		else:
			prepareEnvironment()
			# Start the GUI.
			main()