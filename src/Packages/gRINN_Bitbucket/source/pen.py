#!/usr/bin/env python3

import networkx as nx 
from prody import *
import numpy as np
import pyprind
import argparse
import pandas
import os

def getKongKarplusNetwork(resCorrFile,pdb,resMeanIntEnFile=False,includeCovalents=True,
	corrCutoff=0,intEnCutoff=0,outName=False):

	# Get the number of residues
	sys = parsePDB(pdb)
	numResidues = sys.numResidues()

	# Create the network, and add nodes representing all residues in the protein
	network = nx.Graph()

	for i in range(0,numResidues):
		network.add_node(int(i+1))

	# Load the resCorrFile and determine the maximum edge weight
	resCorrMat = np.loadtxt(resCorrFile)
	resCorrArray = np.squeeze(resCorrMat)
	maxResCorr = np.max(resCorrArray)

	# Load average residue interaction matrix, if provided.
	if resMeanIntEnFile:
		resIntEnMat = np.loadtxt(resMeanIntEnFile)
		if len(resIntEnMat) != numResidues:
			print('Critical error: number of residues in PDB file does not match the size of mean interaction \
				energy matrix')

	# If covalent bonds are also requested, then
	# Connect covalently bound residues with edge distance of -log(maxResCorr)
	if includeCovalents:
		for i in range(0,numResidues-1):
			res1 = sys.select('resindex %i' % i)
			res2 = sys.select('resindex %i' % (i+1))
		
			# Only connect if the two residues are in the same chain
			if (res1.getChids()[0] == res2.getChids()[0]) and (res1.getSegindices()[0] == res2.getSegindices()[1]):
				network.add_edge(int(i+1),int(i+2),weight=float(maxResCorr),distance=0)

	# Add edges between significant interactions. (above intEnCutoff)
	# Assign weights from resCorrMat!
	if resMeanIntEnFile:
		for i in range(0,numResidues):
			for j in range(0,numResidues):
				if not network.has_edge(i+1,j+1):
					if np.abs(resIntEnMat[i][j]) > intEnCutoff:
						network.add_edge(i+1,j+1,weight=float(resCorrMat[i,j]),distance=float(maxResCorr)-float(resCorrMat[i,j]))

	if outName:
		# Write the network to several file formats readable by network analysis packages?
		nx.write_gml(network,outName+'KongKarplus'+'.gml')

	return network

def getRibeiroOrtizNetwork(pdb,resMeanIntEnFile=False,includeCovalents=True,intEnCutoff=1,rmsdEn=5,
	outName='resNetwork'):

	# Get the number of residues
	sys = parsePDB(pdb)
	numResidues = sys.numResidues()

	# Create the network, and add nodes representing all residues in the protein
	network = nx.Graph()

	for i in range(0,numResidues):
		network.add_node(i+1)

	# Load average residue interaction matrix.
	resIntEnMat = np.loadtxt(resMeanIntEnFile)

	# Determine RMSD of interaction energies (needed later on)
	resIntEnArray = [i for i in np.reshape(resIntEnMat,(1,numResidues**2))[0]]
	#rmsdIntEn = np.sqrt(np.mean((resIntEnArray - np.mean(resIntEnArray)) ** 2))

	# Construct an matrix to make edge weights later on according to Ribeiro et al. (2014)
	#X = 0.5*(1-(resIntEnMat-np.mean(resIntEnArray))/(5*rmsdIntEn))

	# Construct an edge weights matrix.
	resIntEnMatNegFavor = np.zeros(np.shape(resIntEnMat))
	for i in range(0,np.shape(resIntEnMat)[0]):
		for j in range(0,np.shape(resIntEnMat)[0]):
			if resIntEnMat[i,j] < 0:
				resIntEnMatNegFavor[i,j] = np.abs(resIntEnMat[i,j])
			else:
				resIntEnMatNegFavor[i,j] = 0

	X = np.abs(resIntEnMatNegFavor)/np.max(np.abs(resIntEnMatNegFavor))

	for i in range(0,numResidues):
		for j in range(0,numResidues):

			if X[i,j] > 0.99:
				X[i,j] = 0.99

			if X[i,j] < 0:
				X[i,j] = 0
				print('alert: negative weight, reassigning to zero.')

	# If covalent bonds are also requested, then
	# Connect covalently bound residues with edge weight of 0.99 and an edge distance of 1/weight.
	if includeCovalents:
		for i in range(0,numResidues-1):
			res1 = sys.select('resindex %i' % i)
			res2 = sys.select('resindex %i' % (i+1))
		
			# Only connect if the two residues are in the same chain
			# Weights as is
			if (res1.getChids()[0] == res2.getChids()[0]) and (res1.getSegindices()[0] == res2.getSegindices()[1]):
				network.add_edge(i+1,i+2,weight=X[i,i+1],distance=1-float(X[i,i+1]))


			# Weights with -log
			#if (res1.getChids()[0] == res2.getChids()[0]) and (res1.getSegindices()[0] == res2.getSegindices()[1]):
			#	network.add_edge(i+1,i+2,{'distance':-np.log(0.99)})

	progbar = pyprind.ProgBar(385**2)

	for i in range(0,numResidues):
		for j in range(0,numResidues):

			if not includeCovalents:
			# Check again for covalent connection. If we are iterating over a residue pair covalently connected, then skip
				if abs(i-j) == 1:
					continue

		# Check whether edges exist between these residues.
			if not network.has_edge(i+1,j+1):

				# Check whether the mean interaction energy between the two residues is above the cutoff value.
				# If yes, continue.
				if abs(float(resIntEnMat[i,j])) >= float(abs(intEnCutoff)):

					# Add an edge between the two residues. Specify distance according to Ribeiro et al. (2014)
					# Also, consider edge weights lower than 0.01 disconnected
					# (again, Ribeiro et al. 2014)
					if X[i,j] < 0.01:
						continue

					# Connect the two residues with edge weight as calculated above and an edge distance of 1/weight.
					# weights as is
					network.add_edge(i+1,j+1,weight=X[i,j],distance=1-float(X[i,j]))
					# weights as -log
					#network.add_edge(i+1,j+1,{'distance':-np.log(X[i,j])})

		progbar.update()

	if outName:
		# Write the network to several file formats readable by network analysis packages?
		nx.write_gml(network,os.path.join(outName,'network'+'.gml'))

	return network

def getProEnNet(inFolder=False,resMeanIntEnFile=False,resCorrFile=False,includeCovalents=True,intEnCutoff=1,
	resCorrCutoff=0.4,outPrefix=False):
	
	if inFolder:
		resMeanIntEnFile = os.path.join(inFolder,'energies_intEnMeanTotal.dat')
		pdb = os.path.join(inFolder,'system_dry.pdb')
		if os.path.exists(os.path.join(inFolder,'energies_resCorr.dat')):
			resCorrFile = os.path.join(inFolder,'energies_resCorr.dat')

	if not resMeanIntEnFile:
		print('You have to provide the mean interaction energy file.')
		raise SystemExit(0)

	networkRO = getRibeiroOrtizNetwork(pdb=pdb,resMeanIntEnFile=resMeanIntEnFile,includeCovalents=includeCovalents,
		intEnCutoff=intEnCutoff,outName=outPrefix)

	networkKK = False
	if not resCorrFile:
		pass
	else:
		pass
		# UNVALIDATED NETWORK ACTIVATION EXPERIMENTAL.
		# getKongKarplusNetwork(resCorrFile=resCorrFile,pdb=pdb,
		# 	resMeanIntEnFile=resMeanIntEnFile,includeCovalents=includeCovalents,corrCutoff=resCorrCutoff,
		# 	intEnCutoff=intEnCutoff,outName=outPrefix)

	return networkRO, networkKK

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
		description='Construct Protein Energy Networks from getResIntEn.py output')

	# Overriding convert_arg_line_to_args in the input parser with our own function.
	parser.convert_arg_line_to_args = convert_arg_line_to_args

	parser.add_argument('--pdb',type=str,nargs=1,help='Path to the PDB file of the protein system')

	parser.add_argument('--resmeanintenfile',type=str,default=[False],nargs=1,
		help='Path to the average interaction energy matrix produced by getResIntEn.py')

	parser.add_argument('--rescorrfile',type=str,default=[False],nargs=1,
		help='Residue correlation matrix produced by getResIntEn.py')

	parser.add_argument('--includecovalents',action='store_true',default=False,
		help='Whether to include covalent bonds or not while constructing the network edges.')

	parser.add_argument('--intencutoff',type=float,nargs=1,default=[1],
		help='Mean (average) interaction energy cutoff when constructing the Ribeiro-Ortiz network \
		(kcal/mol). If an interaction energy time series absolute average value is below this \
		cutoff, not edge will be added between the two residues. By default, the cutoff is 1 kcal/mol.')

	parser.add_argument('--rescorrcutoff',type=float,nargs=1,default=[0.4],
		help='Residue correlation cutoff for inserting the edges between two residues in Kong-Karplus \
		network.')

	parser.add_argument('--infolder',type=str,nargs=1,default=[False],
		help='Output folder specified when calling getResIntEn.py')

	parser.add_argument('--outprefix',type=str,nargs=1,default=['energies_resNetwork_'],
		help='Output file name prefix')

	# Parse input arguments
	args = parser.parse_args()

	pdb = args.pdb[0]

	resMeanIntEnFile = args.resmeanintenfile[0]
	intEnCutoff = float(args.intencutoff[0])
	includeCovalents = args.includecovalents
	resCorrCutoff = float(args.rescorrcutoff[0])
	resCorrFile = args.rescorrfile[0]
	outPrefix = args.outprefix[0]

	inFolder = args.infolder[0]

	_,_ = getProEnNet(inFolder=inFolder,resMeanIntEnFile=resMeanIntEnFile,resCorrFile=resCorrFile,
		includeCovalents=includeCovalents,intEnCutoff=intEnCutoff,resCorrCutoff=resCorrCutoff,
		outPrefix=outPrefix)