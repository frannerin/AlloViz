#!/usr/bin/env python
import re, os, pexpect, panedr, pandas, time, sys, os, psutil
import numpy as np
from prody import *
import logging

class parameters(object):
	def __init__(self):
		f = open(os.path.join(
			os.path.dirname(os.path.abspath(__file__)),'VERSION'),'r')
		self.__version__ = f.readlines()[0]
		f.close()
		self.pdb = None
		self.tpr = None
		self.top = None
		self.traj = None
		self.numCores = None
		self.dielectric = None
		self.sel1 = None
		self.sel2 = None
		self.pairFilterCutoff = None
		self.pairFilterPercentage = None
		self.cutoff = None
		self.stride = None
		self.frameRange = None
		self.namd2exe = None
		self.gmxexe = None
		self.parameterFile = None
		self.calcCorr = None
		self.outFolder = None
		self.dataType = None
		self.logFile = None
		self.logger = None
		self.pairsFiltered = None
		self.pool = None
		self.calcState = 'Idle'

# A method to get the path to the resource for pyinstaller onefile executables.
def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.dirname(os.path.abspath(__file__))

    return os.path.join(base_path, relative_path)

# A method to check whether a file has a handle on it.
def has_handle(fpath):
	for proc in psutil.process_iter():
		try:
			for item in proc.open_files():
				if fpath == item.path:
					return True
		except Exception:
			pass

	return False

# A method to check whether a given command is an executable
def which(program):
	def is_exe(fpath):
		return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

	fpath, fname = os.path.split(program)
	if fpath:
		if is_exe(program):
			return program
	else:
		for path in os.environ["PATH"].split(os.pathsep):
			exe_file = os.path.join(path, program)
			if is_exe(exe_file):
				return exe_file

	return None

def isMemoryEnough(params,traj):
	# Check whether the system has enough memory for multiple processing of the DCD
	trajStats = os.stat(params.traj)
	size = trajStats.st_size

	memory = psutil.virtual_memory()
	if size*params.numCores > memory.available*1.1:
		message = 'System does not have enough memory to handle the computation. '
		'Please either decrease the number of processors (numCores) or increase '
		'the trajectory stride parameter. Aborting now.'
		return False, message
	else:
		return True, "Success"

# A method to get a string containing chain ID, residue name and residue number
# given a ProDy parsed PDB Atom Group and the residue index
def getChainResnameResnum(pdb,resIndex):
	# Get a string for chain+resid+resnum when supplied the residue index.
	selection = pdb.select('resindex %i' % resIndex)
	chain = selection.getChids()[0]
	resName = selection.getResnames()[0]
	resNum = selection.getResnums()[0]
	string = chain+resName+str(resNum)
	return string

# A reverse of the above method
def getResindex(pdb,chainResnameResnum):
	# Get the residue index of a chain resname resnum string.
	matches = re.search('(\D+)(\D{3})(\d+)',chainResnameResnum)
	if matches:
		chain = matches.groups()[0]
		resName = matches.groups()[1]
		resNum = int(matches.groups()[2])
		selection = pdb.select('chain %s and resnum %i' % (chain,resNum))
		resIndex = selection.getResindices()[0]
		return resIndex

# A method that parses the energy log output of NAMD.
def parseEnergiesSingleCoreNAMD(args):

	filePaths = args[0]
	pdb = args[1]
	logFile = args[2]
	logging.basicConfig(format='%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s',
		datefmt='%d-%m-%Y:%H:%M:%S',level=logging.DEBUG,filename=logFile)
	logger = logging.getLogger(__name__)
	# Start a dictionary for storing residue-pair energy values

	energiesDict = dict()
	for filePath in filePaths:
		logger.info('Parsing: '+filePath)
		# Get the interaction residues
		matches = re.search('(\d+)_(\d+)_energies.log',filePath)
		if not matches:
			continue 

		# Get residue indices
		res1 = int(matches.groups()[0])
		res2 = int(matches.groups()[1])

		system = parsePDB(pdb)
		# Get chain-resname-resnum strings
		res1_string = getChainResnameResnum(system,res1)
		res2_string = getChainResnameResnum(system,res2)

		# Read in the first line (header) output file and count number of total lines.
		f = open(filePath,'r')
		lines = f.readlines()

		# Ignore lines not starting with ENERGY:
		lines = [line for line in lines if line.startswith('ENERGY:')]
		f.close()

		lines = [line.strip('\n').split() for line in lines if line.startswith('ENERGY:')]
		lines = [[float(integer) for integer in line[1:]] for line in lines]

		headers = ['Frame','Elec','VdW','Total']
		headerColumns = [0,5,6,10] # Order in which the headers appear in NAMD2 log
		# Frame: 0, Elec: 5, VdW: 6, Total: 10
		numTitles = len(headers)
		# Assign each column into appropriate key's value in energyOutput dict.
		energyOutput = dict()
		for i in range(0,numTitles):
			energyOutput[headers[i]] = [line[headerColumns[i]] for line in lines]

		# Puts this energyOutput dict into energies dict with keys as residue ids
		energiesDict[res1_string+'-'+res2_string] = energyOutput
		# Also store it as res2,res2 (it is the same thing after all)
		energiesDict[res2_string+'-'+res1_string] = energyOutput

	return energiesDict

# A method that parses GMX edr file content (energy output)
def parseEnergiesGMX(gmxExe,pdb,outputFolder,pairsFilteredChunks,edrFiles,logger):

	system = parsePDB(pdb)
	system_dry = system.select('protein or nucleic')
	system_dry = system_dry.select('not resname SOL')
	#gmxExe = '/usr/local/gromacs/bin/gmx' # TEMPORARY!
	# Parse the resulting interact.edr file from the output directory

	# WITHOUT PANEDR
	# proc = pexpect.spawnu('%s energy -f %s/interact.edr -o %s/interact.xvg' % (gmxExe,outputFolder,outputFolder))
	# for pair in pairsFiltered:
	# 	res1 = str(pair[0])
	# 	res2 = str(pair[1])
	# 	proc.send('LJ-SR:res%s-res%s' % (res1,res2))
	# 	proc.sendline()
	# 	proc.send('Coul-SR:res%s-res%s' % (res1,res2))
	# 	proc.sendline()
	# proc.sendline()
	# proc.sendline()
	# proc.kill(1)

	# WITH PANEDR parse edr files to pandas dataframes
	logger.info('Parsing GMX energy output... This may take a while...')
	df = panedr.edr_to_df(os.path.join(outputFolder,'interact0.edr'))
	logger.info('Parsed 1 EDR file.')
	for i in range(1,len(edrFiles)):
		edrFile = edrFiles[i]
		df_pair = panedr.edr_to_df(edrFile)
		# It is possible that an interaction column has been calculated already in accumulated df.
		# This is because of the way GMX energy calculations work.
		# In this case we have to remove them.
		df_pair_columns = df_pair.columns
		df_pair = df_pair[list(set(df_pair_columns)-set(df.columns))]

		df = pandas.concat([df,df_pair],axis=1)
		logger.info('Parsed %i out of %i EDR files...' % (i+1,len(edrFiles)))

	# logger.info('Removing duplicates from energy data frame...')
	# df.to_csv(outputFolder+'/energies.csv')
	# df = pandas.read_csv(outputFolder+'/energies.csv')
	# col2remove = list()
	# for column in df.columns:
	# 	logger.info('Current column: %s' % column)
	# 	matches = re.search('.*\.\d+',column)
	# 	if matches:
	# 		col2remove.append(column)
	# 		logger.info('Will remove column from energy data frame.')

	# df = df.drop(col2remove,1)
	# logger.info('Removing duplicates from energy data frame... Done')
	# os.remove(outputFolder+'/energies.csv')

	# for column in df_pair.columns:
	# 	if column in df.columns:
	# 		df_pair = df_pair.drop(column,axis=1)
	# 		# TEMP
	# 		logger.info('Removed a column that already exists in energies data frame.')

	logger.info('Collecting energy results...')
	energiesDict = dict()
	for i in range(0,len(pairsFilteredChunks)):
		pairsFilteredChunk = pairsFilteredChunks[i]
		energiesDictChunk = dict()
		for pair in pairsFilteredChunk:
			res1_string = getChainResnameResnum(system_dry,pair[0])
			res2_string = getChainResnameResnum(system_dry,pair[1])
			energyDict = dict()

			# Lennard-Jones Short Range interaction
			column_stringLJSR1 = 'LJ-SR:res%i-res%i' % (pair[0],pair[1])
			column_stringLJSR2 = 'LJ-SR:res%i-res%i' % (pair[1],pair[0])
			if column_stringLJSR1 in df.columns:
				column_stringLJSR = column_stringLJSR1
				
			elif column_stringLJSR2 in df.columns:
				column_stringLJSR = column_stringLJSR2

			else:
				logger.exception('At least one required residue interaction was not found in the pair interaction '
					'energy output. Please contact the developer.')
				return

			# Lennard-Jones 1-4 interaction
			column_stringLJ141 = 'LJ-14:res%i-res%i' % (pair[0],pair[1])
			column_stringLJ142 = 'LJ-14:res%i-res%i' % (pair[1],pair[0])
			if column_stringLJ141 in df.columns:
				column_stringLJ14 = column_stringLJ141
				
			elif column_stringLJ142 in df.columns:
				column_stringLJ14 = column_stringLJ142

			else:
				logger.exception('At least one required residue interaction was not found in the pair interaction '
					'energy output. Please contact the developer.')
				return

			# Coulombic Short Range interaction
			column_stringCoulSR1 = 'Coul-SR:res%i-res%i' % (pair[0],pair[1])
			column_stringCoulSR2 = 'Coul-SR:res%i-res%i' % (pair[1],pair[0])
			if column_stringCoulSR1 in df.columns:
				column_stringCoulSR = column_stringCoulSR1
				
			elif column_stringCoulSR2 in df.columns:
				column_stringCoulSR = column_stringCoulSR2

			else:
				logger.exception('At least one required residue interaction was not found in the pair interaction '
					'energy output. Please contact the developer.')
				return

			# Coulombic Short Range interaction
			column_stringCoul141 = 'Coul-14:res%i-res%i' % (pair[0],pair[1])
			column_stringCoul142 = 'Coul-14:res%i-res%i' % (pair[1],pair[0])
			if column_stringCoul141 in df.columns:
				column_stringCoul14 = column_stringCoul141
				
			elif column_stringCoul142 in df.columns:
				column_stringCoul14 = column_stringCoul142

			else:
				logger.exception('At least one required residue interaction was not found in the pair interaction '
					'energy output. Please contact the developer.')
				return

			# Remember that gmx uses the SI units.
			# In terms of energy, this is usually kJ/mol
			# We should convert to kcal/mol to be consistent with NAMD units (and with ourselves)
			kj2kcal = 0.239005736
			enLJSR = np.asarray(df[column_stringLJSR].values)*kj2kcal
			enLJ14 = np.asarray(df[column_stringLJ14].values)*kj2kcal
			enLJ = [enLJSR[j]+enLJ14[j] for j in range(0,len(enLJSR))]
			energyDict['VdW'] = enLJ

			enCoulSR = np.asarray(df[column_stringCoulSR].values)*kj2kcal
			enCoul14 = np.asarray(df[column_stringCoul14].values)*kj2kcal
			enCoul = [enCoulSR[j]+enCoul14[j] for j in range(0,len(enCoulSR))]
			energyDict['Elec'] = enCoul

			energyDict['Total'] = [energyDict['VdW'][j]+energyDict['Elec'][j] for j in range(0,len(energyDict['VdW']))]

			key1 = res1_string+'-'+res2_string
			key1 = key1.replace(' ','')
			key2 = res2_string+'-'+res1_string
			key2 = key2.replace(' ','')
			energiesDictChunk[key1] = energyDict
			energiesDictChunk[key2] = energyDict

		energiesDict.update(energiesDictChunk)
		logger.info('Collected %i out of %i results' % (i+1,len(pairsFilteredChunks)))

	return energiesDict

# A method that creates MDP and NDX file(s) necessary for GMX energy interaction calculation.
def makeNDXMDPforGMX(gmxExe='gmx',pdb=None,tpr=None,soluteDielectric=1,pairsFiltered=None,sourceSel=None,
	targetSel=None,outFolder=os.getcwd(),logger=None):
	
	system = parsePDB(pdb)

	# Modify atom serial numbers to account for possible PDB files with more than 99999 atoms
	system.setSerials(np.arange(1,system.numAtoms()+1))
	
	system_dry = system.select('protein or nucleic')
	system_dry = system_dry.select('not resname SOL')

	if not pairsFiltered and sourceSel and targetSel: # For use outside of getResIntEn

		# Get the source & target selection.
		source = system_dry.select(sourceSel)
		target = system_dry.select(targetSel)

		# For each individual residue included in this selection, get serials, dump them in a dictionary
		# Divide atom serials into chunks of size 15 (this is how index files are made by GMX apparently)
		sourceResIndices = np.unique(source.getResindices())
		sourceResSerials= dict()

		for index in sourceResIndices:
			residue = source.select('resindex %i' % index)
			lenSerials = len(residue.getSerials())
			if lenSerials > 14:
				residueSerials = residue.getSerials()
				sourceResSerials[index] = [residueSerials[i:i+14] for i in range(0,lenSerials,14)]
			else:
				sourceResSerials[index] = np.asarray([residue.getSerials()])


		# Do the same stuff for the target selection.
		# Get the source selection.
		targetResIndices = np.unique(target.getResindices())
		targetResSerials = dict()

		for index in targetResIndices:
			residue = target.select('resindex %i' % index)
			lenSerials = len(residue.getSerials())
			if lenSerials > 14:
				residueSerials = residue.getSerials()
				targetResSerials[index] = [residueSerials[i:i+14] for i in range(0,lenSerials,14)]
			else:
				targetResSerials[index] = np.asarray([residue.getSerials()])

		# Merge the two dicts.
		allSerials = dict()
		allSerials.update(sourceResSerials)
		allSerials.update(targetResSerials)

	else:

		indicesFiltered = np.unique(np.hstack(pairsFiltered))
		allSerials = dict()

		for index in indicesFiltered:
			residue = system_dry.select('resindex %i' % index)
			lenSerials = len(residue.getSerials())
			if lenSerials > 14:
				residueSerials = residue.getSerials()
				allSerials[index] = [residueSerials[i:i+14] for i in range(0,lenSerials,14)]
			else:
				allSerials[index] = np.asarray([residue.getSerials()])

	# Write a standart .ndx file for GMX
	filename = str(outFolder)+'/interact.ndx'

	proc = pexpect.spawnu('%s make_ndx -f %s -o %s' % (gmxExe,tpr,filename))
	proc.logfile = sys.stdout
	proc.expect(u'>.*')
	proc.send('q')
	proc.sendline()
	proc.sendline()

	time.sleep(5)# proc.wait() is a better option here, but it apparently does not work
	# out of the box in Mac OSX. so using time.sleep(2) here.
	proc.kill(1)

	# Append our residue groups to this standart file!
	f = open(filename,'a')
	for key in allSerials:
		f.write('[ res%i ]\n' % key)
		if type(allSerials[key][0]).__name__ == 'ndarray':
			for line in allSerials[key][0:]:
				f.write(' '.join(list(map(str,line)))+'\n')
		else:
			f.write(' '.join(list(map(str,allSerials)))+'\n')
	#f.close()

	# Write the .mdp files necessary for GMX
	mdpFiles = list()

	# Divide pairsFiltered into chunks so that each chunk does not contain
	# more than 200 unique residue indices. This is because GMX does not 
	# allow more than 254 energy groups to be defined in an MDP file and
	# each residue has to be an energy group in our case.
	pairsFilteredChunks = list()
	if len(np.unique(np.hstack(pairsFiltered))) <= 60:
		pairsFilteredChunks.append(pairsFiltered)
	else:
		i = 2
		maxNumRes = len(np.unique(np.hstack(pairsFiltered)))
		while maxNumRes >= 60:
			pairsFilteredChunks = np.array_split(pairsFiltered,i)
			chunkNumResList = list()
			for chunk in pairsFilteredChunks:
				chunkNumResList.append(len(np.unique(np.hstack(chunk))))

			maxNumRes = np.max(chunkNumResList)
			i += 1

	# TEMPORARY
	# Check whether pairsFilteredChunks missed any residue included in pairsFiltered
	
	for pair in pairsFiltered:
		if pair not in np.vstack(pairsFilteredChunks):
			logger.exception('Missing at least one residue in filtered residue pairs. Please contact the developer.')

	# END OF TEMPORARY SECTION

	i = 0
	for chunk in pairsFilteredChunks:
		filename = str(outFolder)+'/interact'+str(i)+'.mdp'
		f = open(filename,'w')
		#f.write('cutoff-scheme = group\n')
		f.write('cutoff-scheme = Verlet\n')
		#f.write('epsilon-r = %f\n' % soluteDielectric)

		chunkResidues = np.unique(np.hstack(chunk))

		resString = ''
		for res in chunkResidues:
			resString += 'res'+str(res)+' '

		#resString += ' SOL'

		f.write('energygrps = '+resString+'\n')

		# Add energygroup exclusions.
		#energygrpExclString = 'energygrp-excl ='

		# GOTTA COMMENT OUT THE FOLLOWING DUE TO TOO LONG LINE ERROR IN GROMPP
		# for key in allSerials:
		# 	energygrpExclString += ' res%i res%i' % (key,key)

		#energygrpExclString += ' SOL SOL'
		#f.write(energygrpExclString)

		f.close()
		mdpFiles.append(filename)
		i += 1

	return mdpFiles, pairsFilteredChunks