#!/usr/bin/env python
from prody import *
import numpy as np
import mdtraj, multiprocessing, pexpect, sys, itertools, argparse, os, pyprind, subprocess, \
re, pickle, types, logging, datetime, psutil, signal, time, pandas, glob, platform, \
traceback, click, copy
from shutil import copyfile, rmtree
from .common import * # Changed by frannerin
from . import corr # Changed by frannerin

# Method for calculating mean interaction energies.
def getResIntEnMean(intEnPickle,pdb,frameRange=False,prefix=''):

	# Load interaction energy pickle file
	intEnFile = open(intEnPickle,'rb')
	intEn = pickle.load(intEnFile)
	numFrames = len(intEn[list(intEn.keys())[0]]['Total'])

	if not frameRange:
		frameRange = [0,numFrames]

	# Get number of residues
	system = parsePDB(pdb)
	system_dry = system.select('protein or nucleic')
	system_dry = system_dry.select('not resname SOL')
	numResidues = len(np.unique(system_dry.getResindices()))

	# Start interaction energy variables
	intEnDict = dict()
	intEnDict['Elec'] = np.zeros((numResidues,numResidues))
	intEnDict['Frame'] = np.zeros((numResidues,numResidues))
	intEnDict['Total'] = np.zeros((numResidues,numResidues))
	intEnDict['VdW'] = np.zeros((numResidues,numResidues))

	progbar = pyprind.ProgBar(numResidues*numResidues)

	filteredButNoInt = list() # Accumulate interactions that were included in calculation but resulted in zero interaction energy.
	for i in range(numResidues):
		i_chainResnameResnum = getChainResnameResnum(system_dry,i)
		for j in range(numResidues):
			j_chainResnameResnum = getChainResnameResnum(system_dry,j)
			keyString = i_chainResnameResnum+'-'+j_chainResnameResnum
			if keyString in intEn:
				intEnDict['Elec'][i,j] = np.mean(intEn[keyString]['Elec'][frameRange[0]:frameRange[1]])
				intEnDict['Elec'][j,i] = np.mean(intEn[keyString]['Elec'][frameRange[0]:frameRange[1]])
				totalMeanEn = np.mean(intEn[keyString]['Total'][frameRange[0]:frameRange[1]])
				intEnDict['Total'][i,j] = totalMeanEn
				intEnDict['Total'][j,i] = totalMeanEn
				intEnDict['VdW'][i,j] = np.mean(intEn[keyString]['VdW'][frameRange[0]:frameRange[1]])
				intEnDict['VdW'][j,i] = np.mean(intEn[keyString]['VdW'][frameRange[0]:frameRange[1]])

				if not totalMeanEn:
					filteredButNoInt.append(keyString)

			else:
				intEnDict['Elec'][i,j] = int(0)
				intEnDict['Elec'][j,i] = int(0)
				intEnDict['Total'][i,j] = int(0)
				intEnDict['Total'][j,i] = int(0)
				intEnDict['VdW'][i,j] = int(0)
				intEnDict['VdW'][j,i] = int(0)

			progbar.update()

	# Save to text
	np.savetxt('%s_intEnMeanTotal.dat' % prefix,intEnDict['Total'])
	np.savetxt('%s_intEnMeanVdW.dat' % prefix,intEnDict['VdW'])
	np.savetxt('%s_intEnMeanElec.dat' % prefix,intEnDict['Elec'])

	# Save in column format as well (only Totals for now)
	f = open('%s_intEnMeanTotal' % prefix+'List.dat','w')
	for i in range(0,len(intEnDict['Total'])):
		for j in range(0,len(intEnDict['Total'][i])):
			value = intEnDict['Total'][i,j]
			if value: # i.e. if it's not equal to zero (included in filtering step or included but was zero)
				f.write('%s\t%s\t%s\n' % (getChainResnameResnum(system_dry,i),getChainResnameResnum(system_dry,j),str(value)))

	f.close()
	
	return intEnDict, filteredButNoInt

# Method for checking whether structure is dry (NAMD-type input)
def isStructureDry(pdb,psf):
	# Load the PDB and PSF files
	pdb = parsePDB(pdb)
	pdbProtein = pdb.select('protein')
	pdbNonProtein = pdb.select('not protein')

	psf = parsePSF(psf)
	psfNonProtein = psf.select('not protein')

	if not pdbNonProtein == None and not psfNonProtein == None:
		return False
	else:
		return True

# Method for preparing input files for NAMD-type data.
def prepareFilesNAMD(params):
	# Detect whether there are non-protein components in the system.
	params.logger.info('Checking whether the structure has non-protein atoms...')
	dryStructure = isStructureDry(params.pdb,params.top)

	if not dryStructure:
		params.logger.info('Non-protein atoms detected in structure...')
		params.logger.info('There are non-protein elements in your input files. Please '
			'consider generating PSF/PDB/DCD files containing only the protein in your '
			'structure.') # Maybe add here a link where the user can check his/her options.

		# Changed by frannerin
		#if sys.stdin.isatty():
		#	if not click.confirm('Do you want to continue?', default=True):
		#		errorSuicide(params,'User requested abort. Aborting now.',removeOutput=False)

	# Proceeding.
	# Just copy psf and pdb files and the trajectory with stride to output folder.
	copyfile(params.pdb,os.path.join(params.outFolder,'system_dry.pdb'))
	copyfile(params.top,os.path.join(params.outFolder,'system_dry.psf'))
	# Load the DCD file, get rid of non-protein sections.
	traj = Trajectory(params.traj)
	pdb = parsePDB(params.pdb)
	traj.link(pdb)
	traj.setAtoms(pdb)
	writeDCD(os.path.join(params.outFolder,'traj_dry.dcd'),
		traj,step=params.stride)
	# Load it back, superpose, save again.
	traj = parseDCD(os.path.join(params.outFolder,'traj_dry.dcd'))
	traj.setAtoms(pdb)
	traj.superpose()
	writeDCD(os.path.join(params.outFolder,'traj_dry.dcd'),traj)

	# Check whether system has enough memory to handle the computation...
	proceed, message = isMemoryEnough(params,os.path.join(params.outFolder,'traj_dry.dcd'))
	if not proceed:
		errorSuicide(params,message,removeOutput=False)

# Method for preparing input files for gromacs-type data.
def prepareFilesGMX(params):
	params.logger.info('Converting TPR to PDB...')

	# Convert tpr to pdb, full system.
	isPDB,messageOut = tpr2pdb(params,params.tpr,
		os.path.join(params.outFolder,'system.pdb'))

	# Check whether file has been created. If not, wait.
	while not os.path.exists(os.path.join(
		params.outFolder,'system.pdb')):
		time.sleep(1)

	# Check whether the file is still being written to...
	while has_handle(os.path.join(
		params.outFolder,'system.pdb')):
		time.sleep(1)

	params.pdb = os.path.join(params.outFolder,'system.pdb')

	params.logger.info('Converting TPR to PDB... Done.')
	pdb = os.path.join(params.outFolder,'system.pdb')
	copyfile(params.tpr,os.path.join(params.outFolder,'system.tpr'))
	tpr = os.path.join(params.outFolder,'system.tpr')

	# Make dry PDB out of the resulting PDB.
	pdb = parsePDB(params.pdb)
	pdbProtein = pdb.select('protein')
	writePDB(os.path.join(params.outFolder,'system_dry.pdb'),pdbProtein)

	# Convert XTC/TRR trajectories to DCD for ProDy compatible analysis...
	params.logger.info('Converting XTC/TRR to DCD...')
	try:
		if params.traj.endswith('.xtc'):
			traj = mdtraj.load_xtc(params.traj,
				top=os.path.join(params.outFolder,'system.pdb'),
				stride=params.stride)
		elif params.traj.endswith('.trr'):
			traj = mdtraj.load_trr(params.traj,
				top=os.path.join(params.outFolder,'system.pdb'),
				stride=params.stride)

		traj.save_trr(os.path.join(params.outFolder,'traj.trr'))
		params.traj = os.path.join(params.outFolder,'traj.trr')

		dataType = 'GMX' # Specify a data type to use later on!
	except:
		params.logger.exception('Could not load the trajectory file provided. Please check your trajectory.')
		return

	traj.save_dcd(os.path.join(params.outFolder,'traj.dcd'))
	# Load back this DCD and continue with it (for code compatibility with ProDy)
	traj = Trajectory(os.path.join(str(params.outFolder),'traj.dcd'))
	traj.link(pdb)
	traj.setAtoms(pdbProtein)

	# Write
	writeDCD(os.path.join(params.outFolder,'traj_dry.dcd'),traj)

	# Load it back and superpose, then write back.
	traj = parseDCD(os.path.join(params.outFolder,'traj_dry.dcd'))
	traj.setAtoms(pdbProtein)
	traj.superpose()
	writeDCD(os.path.join(params.outFolder,'traj_dry.dcd'),traj)
	
	os.remove(os.path.join(params.outFolder,'traj.dcd'))
	params.logger.info('Converting to XTC/TRR to DCD... Done.')

	return params

# Method for calculation of interaction energies on NAMD-type data
# (called by calcEnergiesNAMD)
def calcEnergiesSingleCoreNAMD(args):
	# Input arguments
	pairsFilteredSingleCore = args[0]
	params = args[1]
	psfFilePath = os.path.join(params.outFolder,'system_dry.psf')
	pdbFilePath = os.path.join(params.outFolder,'system_dry.pdb')
	dcdFilePath = os.path.join(params.outFolder,'traj_dry.dcd')
	skip = 1 # We implemented this stride (skip) in the DCD file already.
	pairFilterCutoff = params.pairFilterCutoff
	cutoff = params.cutoff
	environment = 'vacuum'
	soluteDielectric = params.dielectric
	solventDielectric = 80

	outputFolder = os.path.abspath(params.outFolder)
	namd2exe = params.exe
	# paramFile is a list by default, so we should map to get abspath
	paramFile = params.parameterFile
	logFile = os.path.abspath(params.logFile)

	loggingFormat = '%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s'
	logging.basicConfig(format='%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s',
		datefmt='%d-%m-%Y:%H:%M:%S',level=logging.DEBUG,filename=logFile)
	logger = logging.getLogger(__name__)

	# Also print messages to the terminal
	console = logging.StreamHandler()
	console.setLevel(logging.INFO)
	console.setFormatter(logging.Formatter(loggingFormat))
	logger.addHandler(console)

	logger.info('Started an energy calculation thread.')

	# Defining a method to calculate energies in chunks (to show the progress on the screen).
	def calcEnergiesSingleChunk(pairsFilteredSingleChunk,psfFilePath,pdbFilePath,dcdFilePath,skip,
		pairFilterCutoff,cutoff,environment,soluteDielectric,solventDielectric,outputFolder,namd2exe,paramFile,
		logger):

		for pair in pairsFilteredSingleChunk:
			# Write PDB files for pairInteractionGroup specification
			system = parsePDB(pdbFilePath)
			sel1 = system.select(str('resindex %i' % int(pair[0])))
			sel2 = system.select(str('resindex %i' % int(pair[1])))
			# Changing the values of B-factor columns so that they can be recognized by
			# pairInteractionGroup1 parameter in NAMD configuration file.
			sel1.setBetas([1]*sel1.numAtoms())
			sel2.setBetas([2]*sel2.numAtoms())
			pairIntPDB = '%s/%i_%i-temp.pdb' % (outputFolder,pair[0],pair[1])
			pairIntPDB = os.path.abspath(pairIntPDB)
			writePDB(pairIntPDB,system)

			# SAVING ON THE TWO RESIDUE PAIR TO DO LATER ON(NEEDS TESTING)
			#traj = Trajectory(dcdFilePath)
			#traj.link(system)
			
			#traj.setAtoms(system.select('resindex %i %i' % (pair[0],pair[1])))
			#writeDCD('%i_%i-temp.dcd' % (pair[0],pair[1]),traj)
			
			namdConf = '%s/%s_%s-temp.namd' % (outputFolder,pair[0],pair[1])
			f = open(namdConf,'w')

			f.write('structure %s\n' % psfFilePath)
			f.write('paraTypeCharmm on\n')
			if paramFile:
				for file in paramFile:
					#raise SystemExit(0)
					f.write('parameters %s\n' % file)
			else:
				f.write('parameters %s\n' % (sys.path[0]+'/par_all27_prot_lipid_na.inp'))
			f.write('numsteps 1\n')
			f.write('switching off\n')
			f.write('exclude scaled1-4\n')
			f.write('outputname %i_%i-temp\n' % (pair[0],pair[1]))
			f.write('temperature 0\n')
			f.write('COMmotion yes\n')
			f.write('cutoff %d\n' % cutoff)
			
			if environment == 'implicit-solvent':
				f.write('GBIS on\n')
				f.write('solventDielectric %d\n' % solventDielectric)
				f.write('dielectric %d\n' % soluteDielectric)
				f.write('alphaCutoff %d\n' % (float(cutoff)-3)) # Setting GB radius to cutoff for now. We might want to change this behaviour later.
				f.write('SASA on\n')
			elif environment == 'vacuum':
				f.write('dielectric %d\n' % soluteDielectric)
			else:
				f.write('#environment is %s\n' % str(environment))

			f.write('switchdist 10.0\n')
			f.write('pairInteraction on\n')
			f.write('pairInteractionGroup1 1\n')
			f.write('pairInteractionFile %s\n' % pairIntPDB)
			f.write('pairInteractionGroup2 2\n')
			f.write('coordinates %s\n' % pairIntPDB)
			f.write('set ts 0\n')
			#f.write('coorfile open dcd %i_%i-temp.dcd\n' % (pair[0],pair[1]))
			f.write('coorfile open dcd %s\n' % dcdFilePath)
			f.write('while { ![coorfile read] } {\n')
			f.write('\tfirstTimeStep $ts\n')
			f.write('\trun 0\n')
			f.write('\tincr ts 1\n')
		 	#f.write('\tcoorfile skip\n') # Don't need it once you apply stride to tray_dry.dcd in outputfolder
			f.write('}\n')
			f.write('coorfile close')
			f.close()

			# Run namd2 to compute the energies
			try:
				pid_namd2 = subprocess.Popen([namd2exe,namdConf],
					stdout=open(
						os.path.join(outputFolder,'%i_%i_energies.log' % 
							(pair[0],pair[1])),'w'),stderr=subprocess.PIPE)
				_,error = pid_namd2.communicate()
			except KeyboardInterrupt:
				print('Keyboard interrupt detected.')
				sys.exit(0)

			if error:
				#logger.exception('Error while calling NAMD executable:\n'+error).
				error = error.split('\n')
				fatalErrorLine = None

				for i in range(0,len(error)):
					if 'FATAL ERROR:' in error[i]:
						fatalErrorLine = error[i]
						continue

				if fatalErrorLine:
					return fatalErrorLine

			pid_namd2.wait()

			logger.info('Energies saved to %i_%i_energies.log' % (pair[0],pair[1]))
			if not os.path.exists(os.path.join(params.outFolder,'%i_%i_energies.log' % (pair[0],pair[1]))):
				return "gRINN was supposed to generate %i_%i_energies.log but apparently it failed." % (pair[0],pair[1])

		return None

			#subprocess.call('rm %s' % namdConf,shell=True)
			#subprocess.call('rm %s' % pairIntPDB,shell=True)
			#subprocess.call('rm %i_%i-temp*' % (pair[0],pair[1]),shell=True)
			#raise SystemExit(0)
		# Parse the log file and extract necessary energy values

		# Done.

	# Split the pairsFiltered into chunks to print the progress on the screen.
	if len(pairsFilteredSingleCore) >= 100:
		numChunks = 100
	elif len(pairsFilteredSingleCore) < 10:
		numChunks = 1
	else:
		numChunks = 10

	pairsFilteredChunksSingleCore = np.array_split(pairsFilteredSingleCore,numChunks)

	progBar = pyprind.ProgBar(numChunks)

	# Perform the calculations in chunks
	percent = 0

	for pairsFilteredChunk in pairsFilteredChunksSingleCore:
		try:
			errorMessage = calcEnergiesSingleChunk(pairsFilteredChunk,psfFilePath,pdbFilePath,dcdFilePath,skip,
				pairFilterCutoff,cutoff,environment,soluteDielectric,solventDielectric,outputFolder,namd2exe,paramFile,logger)
		except (SystemExit):
			#logger.exception('Fatal error while calling NAMD executable.
			return 'SystemExit'

		if errorMessage:
			return errorMessage

		progBar.update()
		percent = percent + 100/float(numChunks)
		logger.info('Completed calculation percentage: %s' % percent)

	logger.info('Completed a pairwise energy calculation thread.')

# Main method performing interaction energy calculations on NAMD-type data.
def calcEnergiesNAMD(params):
	# Start energy calculation in chunks
	params.logger.info('Splitting the pairs into chunks...')
	params.pairsFilteredChunks = np.array_split(np.asarray(params.pairsFiltered),params.numCores)

	# Define a worker initializer for graceful exit upon ctrl+c
	parent_id = os.getpid()
	def worker_init():
		def sig_int(signal_num, frame):
			print('signal: %s' % signal_num)
			parent = psutil.Process(parent_id)
			for child in parent.children():
				if child.pid != os.getpid():
					print("killing child: %s" % child.pid)
					child.kill()
			print("killing parent: %s" % parent_id)
			parent.kill()
			print("suicide: %s" % os.getpid())
			psutil.Process(os.getpid()).kill()
			os._exit(0)
		signal.signal(signal.SIGINT, sig_int)

	# Catching CTRL+C SIGINT signals.
	def sigint_handler(signum, frame):
		params.logger = logger
		global pool
		pool.terminate()
		pool.join()
		pool.close()
		print('signal: %s' % signum)
		parent = psutil.Process(os.getpid())
		children = parent.children()
		# Loop over children.
		for child in children:
			while child.children():
				# If child has children, note them
				for grandchild in child.children():
					print("killing grandchild: %s" % grandchild.pid)
					try:
						grandchild.kill()
					except:
						pass

			if child.pid != os.getpid():
				print("killing child: %s" % child.pid)
				child.kill()

		#time.sleep(5)
		if sys.stdin.isatty():
			if not click.confirm('Would you like to delete the output folder?', default=True):
				errorSuicide(params,'Keyboard interrupt detected. Aborting now.',removeOutput=False)
			else:
				errorSuicide(params,'Keyboard interrupt detected. Aborting now.',removeOutput=True)
		else:
			errorSuicide(params,'GUI interrupt detected. Aborting now.',removeOutput=False)
		
		#print("killing parent: %s" % parent_id)
		#parent.kill()
		#print("suicide: %s" % os.getpid())
		#psutil.Process(os.getpid()).kill()
		os._exit(0)
	signal.signal(signal.SIGINT, sigint_handler)

	global pool
	pool = multiprocessing.Pool(params.numCores)
	logger = params.logger

	# Use map_aysnc on the previously created multiprocessing pool to spawn multiple singe core
	# energy calculation threads.
	# get(9999999) below is necessary to let the map respond without blocking the spawned threads.
	# This is a python bug in 2.7
	params.logger.info('Starting threads for interaction energy calculation...')
	# Strip logger away from params temporarily to be able to map.
	params.logger = None
	results = pool.map_async(calcEnergiesSingleCoreNAMD,
		zip(params.pairsFilteredChunks,itertools.repeat(params))).get(9999999)
	params.logger = logger
	
	# If the pool return at least one 'SystemExit' string
	# Abort
	# see the calcEnergiesSingleCoreNAMD method)
	
	if 'SystemExit' in results:
		removeOutput = False if sys.stdin.isatty() else False
		errorSuicide(params,'Critical error while calling NAMD executable. \n\n'
			'Error could not be identified in detail. Please inspect your input data carefully.\n'
			'If the error persists, contact us. Aborting now.',
			removeOutput=removeOutput)
	elif results[0] is not None:
		if 'FATAL ERROR: ' in results[0]:
			removeOutput = False if sys.stdin.isatty() else False
			errorMessage = results[0] # Cause with multiple CPUs multiple outputs are possible.
			errorSuicide(params,'Fatal error from NAMD: '+
				errorMessage.lstrip('FATAL ERROR:'),removeOutput=removeOutput)

	# Parse the specified outFolder after energy calculation is done.
	outFolderFileList = os.listdir(params.outFolder)

	energiesFilePaths = list()
	for fileName in outFolderFileList:
		if fileName.endswith('energies.log'):
			energiesFilePaths.append(os.path.join(params.outFolder,fileName))

	energiesFilePathsChunks = np.array_split(list(energiesFilePaths),
		params.numCores)

	parsedEnergiesResults = pool.map_async(parseEnergiesSingleCoreNAMD,
		zip(energiesFilePathsChunks,itertools.repeat(os.path.join(
			params.outFolder,'system_dry.pdb')),
			itertools.repeat(params.logFile))).get(9999999)

	parsedEnergies = dict()
	for parsedEnergiesResult in parsedEnergiesResults:
		parsedEnergies.update(parsedEnergiesResult)

	pool.close()
	pool.join()

	params.parsedEnergies = parsedEnergies
	return params

# Main method performing interaction energy calculations on gromacs-type data.
def calcEnergiesGMX(params):

	params.logger.info('Started an energy calculation thread.')

	# Prevent backup making while calculating energies.
	os.environ["GMX_MAXBACKUP"] = "-1"

	# Make an index and MDP file with the pairs filtered.
	#gmxExe = 'gmx'
	mdpFiles,pairsFilteredChunks = makeNDXMDPforGMX(gmxExe=params.exe,
		pdb=params.pdb,tpr=params.tpr,soluteDielectric=params.dielectric,
		pairsFiltered=params.pairsFiltered,outFolder=params.outFolder,
		logger=params.logger)

	# Call gromacs pre-processor (grompp) and make a new TPR file for each pair and calculate energies for each pair.
	i = 0
	edrFiles = list()
	for i in range(0,len(mdpFiles)):
		mdpFile = mdpFiles[i]
		tprFile = mdpFile.rstrip('.mdp')+'.tpr'
		edrFile = mdpFile.rstrip('.mdp')+'.edr'

		args = [params.exe,'grompp','-f',mdpFile,'-n',
			os.path.join(params.outFolder,'interact.ndx'),'-p',params.top,'-c',
			params.tpr,'-o',tprFile,'-maxwarn','20']
		proc = subprocess.Popen(args)
		proc.wait()

		# Catching CTRL+C SIGINT signals.
		def sigint_handler(signum, frame):
			proc.kill()
			if sys.stdin.isatty():
				if not click.confirm('Would you like to delete the output folder?', default=True):
					errorSuicide(params,'Keyboard interrupt detected. Aborting now.',removeOutput=False)
				else:
					errorSuicide(params,'Keyboard interrupt detected. Aborting now.',removeOutput=True)
			else:
				errorSuicide(params,'GUI interrupt detected. Aborting now.',removeOutput=False)

		signal.signal(signal.SIGINT,sigint_handler)

		proc = subprocess.Popen([params.exe,'mdrun','-rerun',params.traj,'-v','-s',tprFile,
			'-e',edrFile,'-nt',str(params.numCores)],stdout=subprocess.PIPE,stderr=subprocess.PIPE)

		error = proc.communicate()[1]
		if error:
			error = error.split('\n') # Splitting into lines to be able to process each line separately.
			fatalErrorLines = None

			# Collect fatal error and subsequent lines.
			for j in range(0,len(error)):
				if 'Fatal error' in error[j]:
					fatalErrorLines = error[j:]
					continue

			if fatalErrorLines:
				fatalError = '\n'.join(fatalErrorLines)
				message = 'Fatal error from gmx:\n\n' + fatalError
				errorSuicide(params,message)
			# else:
			# 	error = '\n'.join(error)
			# 	message = 'Error from gmx:\n\n' + error
			# 	errorSuicide(params,message)
				
		proc.wait()

		edrFiles.append(edrFile)

		params.logger.info('Completed calculation percentage: '+str((i+1)/float(len(mdpFiles))*100))

	return edrFiles, pairsFilteredChunks

# Method for filtering pairs to include in the calculation based on pairwise distances using
# criteria specified by the user.
def filterPairs(params):
	
	system = parsePDB(os.path.join(params.outFolder,'system_dry.pdb'))
	traj = parseDCD(os.path.join(params.outFolder,'traj_dry.dcd'))

	try:
		sourceCA = system.select(str(params.sel1)+' and name CA')
	except:
		params.logger.exception('Could not select Selection 1 residue group. Aborting now.')
		return

	numSource = len(sourceCA)
	sourceResids = sourceCA.getResindices()
	sourceResnums = sourceCA.getResnums()
	sourceSegnames = sourceCA.getSegnames()

	allResiduesCA = system.select('name CA')
	numResidues = len(allResiduesCA)
	numTarget = numResidues

	# By default, targetResids are all residues.
	targetResids = np.arange(numResidues)
	
	# Get target selection residues
	try:
		targetCA = system.select(str(params.sel2+' and name CA'))
		numTarget = len(targetCA)
		targetResids = targetCA.getResindices()
	except:
		params.logger.exception('Could not select Selection 2 residue group. Aborting now.')
		return

	# Generate all possible unique pairwise residue-residue combinations
	pairProduct = itertools.product(sourceResids,targetResids)
	pairSet = set()
	for x,y in pairProduct:
		if x != y:
			pairSet.add(frozenset((x,y)))

	# Split the pair set list into chunks according to number of cores
	pairChunks = np.array_split(list(pairSet),params.numCores)

	params.logger.info('Starting the filtering step...')

	# Continue with filtering operation
	traj.setAtoms(allResiduesCA)
	coordSets = traj.getCoordsets()

	# Start a contact matrix (Kirchhoff matrix)
	kh = np.zeros((numResidues,numResidues))

	# Accumulate contact matrix as the sim progresses
	calculatedPercentage = 0
	monitor = 0

	for i in range(0,len(coordSets),1):
		coordSet = coordSets[i]
		gnm = GNM('GNM')
		gnm.buildKirchhoff(coordSet,cutoff=params.pairFilterCutoff)
		kh = kh + gnm.getKirchhoff()
		monitor = monitor + 1
		calculatedPercentage = (float(monitor)/float(len(coordSets)))*100
		if calculatedPercentage > 100: calculatedPercentage = 100
		params.logger.info('Filtered pairs percentage: %s' % str(calculatedPercentage))

	# Get whether contacts are below cutoff for the specified percentage of simulation
	pairsInclusionFraction = np.abs(kh)/(len(traj)/float(1))
	pairsFilteredFlag = pairsInclusionFraction > params.pairFilterPercentage*0.01

	pairsFiltered = list()
	#concatSourceTargetResids = np.concatenate([sourceResids,targetResids])
	for sourceResid in sourceResids:
		for targetResid in targetResids:
			if sourceResid == targetResid:
				continue
			elif pairsFilteredFlag[sourceResid,targetResid] > 0:
				pairsFiltered.append(sorted([sourceResid,targetResid]))

	pairsFiltered = sorted(pairsFiltered)
	pairsFiltered = [list(x) for x in set(tuple(x) for x in pairsFiltered)]
	
	# file = open('pairsFiltered.txt','w')
	# for pair in pairsFiltered:
	#  	file.write('%.2f-%i-%i\n' % (pairsInclusionFraction[pair[0],pair[1]],pair[0],pair[1]))
	# file.close()

	if not pairsFiltered:
		params.logger.exception('Filtering step did not yield any pairs. '
			'Either your cutoff value is too small or the percentage value is too high.')
		return

	params.logger.info('Number of interaction pairs selected after filtering step: %i' % len(pairsFiltered))

	params.pairsFiltered = pairsFiltered
	return params

# A method for collecting results from the output directory.
def collectResults(params):
	params.logger.info('Collecting results...')

	# Prepare a pandas data table from parsed energies, write it to new files depending on type of energy
	df_total = pandas.DataFrame()
	df_elec = pandas.DataFrame()
	df_vdw = pandas.DataFrame()
	for key,value in list(params.parsedEnergies.items()):
		df_total[key] = value['Total']
		df_elec[key] = value['Elec']
		df_vdw[key] = value['VdW']

	params.logger.info('Saving results to '+os.path.join(params.outFolder,'energies_intEnTotal.csv'))
	df_total.to_csv(os.path.join(params.outFolder,'energies_intEnTotal.csv'))
	params.logger.info('Saving results to '+os.path.join(params.outFolder,'energies_intEnElec.csv'))
	df_elec.to_csv(os.path.join(params.outFolder,'energies_intEnElec.csv'))
	params.logger.info('Saving results to '+os.path.join(params.outFolder,'energies_intEnVdW.csv'))
	df_vdw.to_csv(os.path.join(params.outFolder,'energies_intEnVdW.csv'))

	params.logger.info('Saving results to '+os.path.join(params.outFolder,'energies.pickle'))
	file = open(os.path.join(params.outFolder,'energies.pickle'),'wb')
	pickle.dump(params.parsedEnergies,file)
	file.close()

	params.logger.info('Getting mean interaction energies...')
	# Save average interaction energies as well!
	intEnDict, filteredButNoInt = getResIntEnMean(os.path.join(params.outFolder,'energies.pickle'),
		os.path.join(params.outFolder,'system_dry.pdb'),
		prefix=os.path.join(params.outFolder,'energies'))

	# Report interactions with zero mean.
	for noint in filteredButNoInt:
		params.logger.info('The interaction %s was included in energy calculation but yielded '
			'zero kcal/mol mean interaction energy.' % noint)	

	return params
	if resIntCorr:
		logger.info('Beginning residue interaction energy correlation calculation...')
		getResIntCorr.getResIntCorr(inFile=os.path.join(
			outputFolder,'energies_intEnTotal.csv'),
			pdb=pdb,meanIntEnCutoff=resIntCorrAverageIntEnCutoff,
			outPrefix=os.path.join(outputFolder,'energies'),logger=logger)

# Cleaning up the output folder.
def cleanUp(params):
	params.logger.info('Cleaning up...')
	# Delete all namd-generated energies file from output folder.
	for item in glob.glob(os.path.join(params.outFolder,'*_energies.log')):
		os.remove(item)

	for item in glob.glob(os.path.join(params.outFolder,'*temp*')):
		os.remove(item)

	# Delete all gromacs-generated energies file from output folder.
	for item in glob.glob(os.path.join(params.outFolder,'interact*')):
		os.remove(item)

	for item in glob.glob(os.path.join(params.outFolder,'*.trr')):
		os.remove(item)

	if os.path.exists(os.path.join(params.outFolder,'traj.dcd')):
		os.remove(os.path.join(params.outFolder,'traj.dcd'))

# A small helper method which removes output directory 
# upon user's request.
def errorSuicide(params,message,removeOutput=False):
	params.logger.exception(message)
	if removeOutput:
		rmtree(params.outFolder,ignore_error=True)
	#psutil.Process(os.getpid()).kill()
	# Exit normally after printing the error to the log file.
	os._exit(0)

# Main method serving as entry point for interaction energy calculations
# on NAMD data.
def calcNAMD(params):
	# Prepare input files for NAMD energy calculation.
	prepareFilesNAMD(params) 

	# Filter pairs.
	params = filterPairs(params)
	pairsFiltered = copy.deepcopy(params.pairsFiltered)

	# Calculate interaction energies.
	params = calcEnergiesNAMD(params)
	if not pairsFiltered == params.pairsFiltered:
		params.logger.error('pairsFiltered do not match...')

	# Collect results.
	params = collectResults(params)

	# Clean up
	cleanUp(params)

	return params

# Main method serving as entry point for interaction energy calculations
# on NAMD data.
def calcGMX(params):
	# Prepare input files for GMX energy calculation.
	params = prepareFilesGMX(params)

	# Filter pairs.
	params = filterPairs(params)

	# Calculate interaction energies.
	edrFiles,pairsFilteredChunks = calcEnergiesGMX(params)

	# Parse resulting energy EDR files.
	params.parsedEnergies = parseEnergiesGMX(gmxExe=params.exe,
		pdb=os.path.join(params.outFolder,'system.pdb'),
		pairsFilteredChunks=pairsFilteredChunks,
		outputFolder=params.outFolder,edrFiles=edrFiles,
		logger=params.logger)

	# Collect results
	params = collectResults(params)

	# Clean up
	cleanUp(params)

	return params

# Method to convert TPR to PDB files.
def tpr2pdb(params,tpr,pdb):
	# Convert tpr to pdb, selecting just Protein.
	# Apparently directly spawning gmx in the following does not work as expect in OSX
	# Prepending bash -c to the command line prior to gmx.
	proc = subprocess.Popen('bash -c "%s editconf -f %s -o %s"' % 
		(params.exe,tpr,pdb),stdout=subprocess.PIPE,
		stderr=subprocess.PIPE,shell=True)
	_,error = proc.communicate()
	#if error:
	#	message = repr(error)
	#	return False, message

	proc.wait()
	start_time = time.time()
	time_elapsed = 0 # Seconds
	while not os.path.exists(pdb) and time_elapsed < 30:
		time.sleep(1) # using time.sleep(X) instead, sleeping for X seconds to let the bg process complete work
		time_elapsed = time.time() - start_time

	if not os.path.exists(pdb):
		message = 'Could not extract PDB out of TPR file. Aborting now.'
		return False, message

	# Check whether the file is still being written to...
	while has_handle(pdb):
		time.sleep(1)

	# Check whether there any chain ids not assigned a valid letter.
	noChid = False
	system = parsePDB(pdb)
	systemProtein = system.select('protein')
	chids = systemProtein.getChids()
	for chid in chids:
		if not chid.strip():
			noChid = True
	
	if noChid:
		if params.logger: # Considering the case when this method is called during argument checking step..
			params.logger.info('At least one atom with no chain IDs present. Assigning the default chain ID P to such atoms right now...')
		system_noChids = system.select('chain " "')
		system_noChids.setChids(['P']*system_noChids.numAtoms())
		writePDB(pdb,system)			

	return True, "Success'"

# Method to check args and get params if they are valid
def getParams(args):

	# Make a new parameters object.
	params = parameters()

	# Check whether the output folder exists. If it exists, abort.
	outFolder = os.path.abspath(args.outfolder[0])
	currentFolder = os.getcwd()
	if outFolder != currentFolder:
		if os.path.exists(outFolder):
			print("The output folder exists. Please delete this folder or "
				" specify a folder path that does not exist. Aborting now.")
			sys.exit(1)
		elif not os.access(os.path.abspath(
			os.path.dirname(outFolder)), os.W_OK):
			print("Can't write to the output folder path. Do you have write access?")
			return
		else:
			params.outFolder = outFolder
			params.logFile = os.path.join(os.path.abspath(outFolder),'grinn.log')

	params.numCores = args.numcores[0]
	frameRange = args.framerange

	if len(frameRange) > 1:
		params.frameRange = np.asarray(frameRange)
	elif len(frameRange) == 1:
		if not frameRange[0]:
			params.frameRange = False
	else:
		message = 'Invalid frame range. Aborting now.'
		return params, False, message

	params.stride = args.stride[0]

	params.dielectric = args.dielectric[0]

	params.pairFilterCutoff = args.pairfiltercutoff[0]

	if params.pairFilterCutoff < 4:
		message = 'Filtering distance cutoff value can not be smaller than 4. Aborting now.'
		return params, False, message

	params.pairFilterPercentage = args.pairfilterpercentage[0]

	params.cutoff = args.cutoff[0]

	if not type(args.sel1) == str:
		if len(args.sel1) > 1:
			params.sel1 = ' '.join(args.sel1)
		else:
			params.sel1 = args.sel1[0]

	if not type(args.sel2) == str:
		if len(args.sel2) > 1:
			params.sel2 = ' '.join(args.sel2)
		else:
			params.sel2= args.sel2[0]

	# Check input simulation data.
	if not args.top[0]:
		message = "You must specify a valid topology file (PSF or TOP). Aborting now."
		return params, False, message
	else:
		params.top = os.path.abspath(args.top[0])
		if params.top.lower().endswith('.psf'):
			try:
				topology = parsePSF(params.top)
			except:
				message = "Could not load your PSF file. Aborting now."
				return params, False, message

	if args.pdb[0] and args.tpr[0]:
		message = "You can't specify a PDB and a TPR file at the same time. Please specify either "
		"a PDB for NAMD data or a TPR for GROMACS data. Aborting now."
		return params, False, message

	if args.pdb[0]:
		try:
			system = parsePDB(os.path.abspath(args.pdb[0]))
			systemProtein = system.select(str('protein or nucleic'))
			params.pdb = os.path.abspath(args.pdb[0])
			params.dataType = 'namd'
		except:
			message = "Could not load your PDB file. Aborting now."
			return params, False, message

		try:
			sysSel1 = system.select(params.sel1)
			sysSel2 = system.select(params.sel2)
		except:
			message = 'Could not select sel1 or sel2 in the PDB file. Aborting now.'
			return params, False, message

		# Check whether there any chain ids not assigned a valid letter.
		chids = systemProtein.getChids()
		for chid in chids:
			if not chid.strip():
				message = 'There is at least one residue with no chain ID assigned to it. This is not '\
				'allowed. Aborting now...'
				return params, False, message

		numResidues = len(np.unique(systemProtein.getResindices()))
		for resindex in np.unique(systemProtein.getResindices()):
			residue = systemProtein.select(str('resindex %i' % resindex))
			index = np.unique(residue.getResnames())
			if len(index) > 1:
				message = 'There are multiple residues with the same residue index in your PDB file. '\
				' This is not allowed. Aborting now...'
				return params, False, message

	elif args.tpr[0]:
		# Unfortunately I don't know of a good way to check valid GMX tpr data.
		if not args.tpr[0].lower().endswith('.tpr'):
			message = "The TPR file must have extension .tpr. Aborting now."
		else:
			params.tpr = os.path.abspath(args.tpr[0])
			params.dataType = 'gmx'

	else:
		message = "Please specify either a PDB for NAMD data or a TPR for GROMACS data. "
		"Aborting now."
		return params, False, message

	if not args.traj[0]:
		message = "You have not specified a trajectory file!"
		return params, False, message
	else:
		params.traj = os.path.abspath(args.traj[0])

	# Check whether given exe is actually an exe!
	# If not, abort.
	if not args.exe[0]:
		message = "You have not specified a NAMD2 or GMX executable!"
		return params, False, message
	if os.path.exists(os.path.join(os.getcwd(),args.exe[0])):
		params.exe = os.path.abspath(args.exe[0])
	else:
		params.exe = args.exe[0]

	isExe = which(params.exe)
	if not isExe:
		message = "NAMD2/GMX exe you specified does not appear to be a valid executable. "
		"Aborting now."
		return params, False, message

	# Check extension combinations.
	_,trajExt = os.path.splitext(params.traj)
	if params.dataType == 'namd':
		_,topExt = os.path.splitext(params.top)
		_,pdbExt = os.path.splitext(params.pdb)
		exts = [topExt.lower(),pdbExt.lower(),trajExt.lower()]
		if exts != ['.psf','.pdb','.dcd']:
			message = 'Invalid PSF/PDB/DCD file extensions. Aborting now.'
			return params, False, message

		try:
			trajectory = Trajectory(params.traj)
		except:
			message = 'Could not load the DCD file provided. Aborting now.'
			return params, False, message

		# Check whether stride is higher than the number of frames in trajectory:
		if params.stride > trajectory.numFrames():
			message = 'Stride value is higher than the number of frames in the trajectory. '\
			'Please use a lower stride value.'
			return params, False, message

		# Check whether a parameter file is supplied.
		parameterFile = args.parameterfile
		for paramFile in parameterFile:
			if not paramFile:
				message = 'You must supply a parameter file for NAMD. Aborting now.'
				return params, False, message

		params.parameterFile = [os.path.abspath(paramFile) for paramFile in parameterFile]
		#print(params.parameterFile)
		#return params, False, "what the hell?"
		
	elif params.dataType == 'gmx':

		if platform.system() == 'Windows':
			message = 'GROMACS data on Windows is not supported. Aborting now.'
			return params, False, message

		_,tprExt = os.path.splitext(params.tpr)
		_,topExt = os.path.splitext(params.top)
		exts = [topExt.lower(),tprExt.lower(),trajExt.lower()]
		if exts != ['.top','.tpr','.trr'] and exts != ['.top','.tpr','.xtc']:
			message = 'Invalid TOP/TPR/XTC/TRR file extensions. Aborting now.'
			return params, False, message

		# Check whether a PDB can be extracted from the TPR.
		isPDB,messageOut = tpr2pdb(params,params.tpr,'dummy.pdb')
		if not isPDB:
			message = 'Could not extract a structure from input TPR.'
			message = message + ' Executable produced the following : ' 
			message = message + messageOut
			return params, False, message
		else:
			try:
				system = parsePDB('dummy.pdb')
				systemProtein = system.select(str('protein or nucleic'))
				os.remove('dummy.pdb')
			except:
				os.remove('dummy.pdb')
				message = 'Could not load the extracted PDB file from TPR. '
				'Aborting now.'
				return params, False, message

			try:
				sysSel1 = system.select(params.sel1)
				sysSel2 = system.select(params.sel2)
			except:
				message = 'Could not select sel1 or sel2 in the PDB file extracted from the '
				'TPR. Aborting now.'
				return params, False, message

	params.calcCorr = args.calccorr

	return params, True, "Success"

# Main method starting the work
def getResIntEn(args):
	
	# Check whether input arguments are valid and get parameters!
	global params
	print('Checking input arguments...')
	params, isArgsValid, message = getParams(args)

	# Create the output folder now so that we can start logging.
	# Creating this file right now is important because the calcGUI 
	# will monitor this file as well.
	try:
		os.makedirs(params.outFolder)
		f = open(params.logFile,'w')
		f.close()
	except:
		print('Failed to create the output directory. Do you have write access?')
		sys.exit(0)

	# Start logging.
	loggingFormat = '%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s'
	logging.basicConfig(format=loggingFormat,datefmt='%d-%m-%Y:%H:%M:%S',level=logging.DEBUG,
		filename=params.logFile)
	params.logger = logging.getLogger(__name__)
	
	# Also print messages to the terminal
	console = logging.StreamHandler()
	console.setLevel(logging.INFO)
	console.setFormatter(logging.Formatter(loggingFormat))
	params.logger.addHandler(console)

	# Check whether the arguments are valid. If not, remove the output folder and abort.
	if not isArgsValid:
		# Check whether the script was called from a terminal.
		if sys.stdin.isatty():
			errorSuicide(params,message,removeOutput=False)
			return
		else:
			errorSuicide(params,message,removeOutput=False)
			return

	params.logger.info('Argument check completed. Proceeding...')

	params.logger.info('Started calculation.')

	# Proceed with the appropriate method depending on the input data type.
	if params.dataType == 'namd':
		params = calcNAMD(params)
	elif params.dataType == 'gmx':
		params = calcGMX(params)

	# Get correlations, if the user requested.
	if params.calcCorr:
		args.corrprefix = [os.path.join(params.outFolder,'energies')]
		args.corrinfile = [os.path.join(params.outFolder,'energies_intEnTotal.csv')]
		args.pdb = [params.pdb]
		corr.getResIntCorr(args,logFile=None,logger=params.logger)

	params.logger.info('FINAL: Computation sucessfully completed. Thank you for using gRINN.')
	return

if __name__ == '__main__':
	print('Please do not call this script directly. Use python grinn.py -calc <arguments> '
		'instead.')
	sys.exit(0)


