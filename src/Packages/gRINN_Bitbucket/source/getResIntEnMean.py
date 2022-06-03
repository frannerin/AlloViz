#!/usr/bin/env python3
from prody import *
import multiprocessing
import numpy as np
import sys, itertools, argparse, os, pyprind, subprocess
import re, pickle, types
from common import getChainResnameResnum

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

	for i in range(numResidues):
		i_chainResnameResnum = getChainResnameResnum(system_dry,i)
		for j in range(numResidues):
			j_chainResnameResnum = getChainResnameResnum(system_dry,j)
			keyString = i_chainResnameResnum+'-'+j_chainResnameResnum
			if keyString in intEn:
				intEnDict['Elec'][i,j] = np.mean(intEn[keyString]['Elec'][frameRange[0]:frameRange[1]])
				intEnDict['Elec'][j,i] = np.mean(intEn[keyString]['Elec'][frameRange[0]:frameRange[1]])
				intEnDict['Total'][i,j] = np.mean(intEn[keyString]['Total'][frameRange[0]:frameRange[1]])
				intEnDict['Total'][j,i] = np.mean(intEn[keyString]['Total'][frameRange[0]:frameRange[1]])
				intEnDict['VdW'][i,j] = np.mean(intEn[keyString]['VdW'][frameRange[0]:frameRange[1]])
				intEnDict['VdW'][j,i] = np.mean(intEn[keyString]['VdW'][frameRange[0]:frameRange[1]])

			else:
				intEnDict['Elec'][i,j] = 0
				intEnDict['Elec'][j,i] = 0
				intEnDict['Total'][i,j] = 0
				intEnDict['Total'][j,i] = 0
				intEnDict['VdW'][i,j] = 0
				intEnDict['VdW'][j,i] = 0

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
			if value: # i.e. it it's not equal to zero
				f.write('%s\t%s\t%s\n' % (getChainResnameResnum(system_dry,i),getChainResnameResnum(system_dry,j),str(value)))

	f.close()
	
	return intEnDict

def convert_arg_line_to_args(arg_line):
	# To override the same method of the ArgumentParser (to read options from a file)
	# Credit and source: hpaulj from StackOverflow
	# http://stackoverflow.com/questions/29111801/using-fromfile-prefix-chars-with-multiple-arguments-nargs#
	for arg in arg_line.split():
		if not arg.strip():
			continue
		yield arg

if __name__ == '__main__':

	# INPUT PARSING
	parser = argparse.ArgumentParser(fromfile_prefix_chars='@',
		description='Calculate the average matrix of previously calculated pairwise nonbonded interaction energies \
		between residues of two selected groups of atoms of a protein over the course of a molecular dynamics DCD trajectory.')

	# Overriding convert_arg_line_to_args in the input parser with our own function.
	parser.convert_arg_line_to_args = convert_arg_line_to_args

	parser.add_argument('--intenpickle',type=str,nargs=1,help='Name of the corresponding pickle energies \
		file')

	parser.add_argument('--pdb',type=str,nargs=1,help='Name of the corresponding pdb file of the system.')

	parser.add_argument('--prefix',type=str,nargs=1,help='Prefix for output files.')

	parser.add_argument('--framerange',type=int,default=[False],nargs='+',help='If specified, then only FRAMERANGE\
		section of the trajectory will be handled')

	# Parsing input arguments
	args = parser.parse_args()

	intenpickle = args.intenpickle[0]
	pdb = args.pdb[0]
	frameRange = args.framerange
	prefix = args.prefix[0]
	
	if len(args.framerange) > 1:
		frameRange = np.asarray(args.framerange)
	else:
		frameRange = args.framerange[0]

	getResIntEnMean(pdb=pdb,intEnPickle=intenpickle,frameRange=frameRange,prefix=prefix)