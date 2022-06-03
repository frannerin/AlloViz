#!/usr/bin/env python
from natsort import natsorted
import numpy as np
import pyprind
import os
import math
import itertools
import multiprocessing
from . import calc # Changed by frannerin
import argparse
import datetime
import logging
import signal
import pandas
import psutil
import re
from prody import *
from .common import getChainResnameResnum # Changed by frannerin

def getResIntCorr(args,logFile=None,logger=None):

	inFile = args.corrinfile[0]
	pdb = args.pdb[0]
	meanIntEnCutoff = args.corrintencutoff[0]
	numCores = args.numcores[0]
	outPrefix = args.corrprefix[0]
	if not logFile:
		logFile = 'grinncorr.log'
	
	if not logger and logFile:
		# Start logging.
		loggingFormat = '%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s'
		logging.basicConfig(format=loggingFormat,datefmt='%d-%m-%Y:%H:%M:%S',level=logging.DEBUG,
			filename=logFile)
		logger = logging.getLogger(__name__)
		
		# Also print messages to the terminal
		console = logging.StreamHandler()
		console.setLevel(logging.INFO)
		console.setFormatter(logging.Formatter(loggingFormat))
		logger.addHandler(console)
	
	logger.info('Started residue interaction energy calculation.')

	logger.info('Reading input CSV...')
	# Read in interaction energy time series from the getResIntEn csv output.
	df = pandas.read_csv(inFile,engine='python')
	logger.info('Reading input CSV... completed.')

	# Get number of residues
	system = parsePDB(pdb)
	systemProtein = system.select('protein or nucleic')
	numResidues = len(np.unique(systemProtein.getResindices()))
	# Convert the interaction energy time series to a 3D matrix
	intEnMat = np.zeros((numResidues,numResidues,len(df)))

	logger.info('Calculating interaction energy correlations...')
	for i in range(0,numResidues):
		for j in range(0,numResidues):
			df_col1 = getChainResnameResnum(system,i)+'-'+getChainResnameResnum(system,j)
			df_col2 = getChainResnameResnum(system,j)+'-'+getChainResnameResnum(system,i)
			if df_col1 in df.columns:
				df_col = df_col1
			elif df_col2 in df.columns:
				df_col = df_col2
			else:
				intEnMat[i,j] = np.zeros(len(df))
				intEnMat[j,i] = np.zeros(len(df))
				continue

			# Only take interactions whose average is above meanIntEnCutoff
			if np.mean(df[df_col].values) >= meanIntEnCutoff:
				intEnMat[i,j] = df[df_col1].values
				intEnMat[j,i] = df[df_col1].values
			else:
				intEnMat[i,j] = np.zeros(len(df))
				intEnMat[j,i] = np.zeros(len(df))

		percentCalculated = ((i+1)/float(numResidues))*25# /4 because this is only quad of calculation.
		logger.info('Interaction energy correlation calculated percentage: %f' % percentCalculated)


	# Calculate pearson product moment correlation between all interaction energy pairs
	# using linear algebra (matrix formalism)
	# This is much more efficient than using for loops

	progbar = pyprind.ProgBar(numResidues)

	# Store correlations in a dictionary.
	sigcorrs = dict()
	for i in range(0,numResidues):
		# First get all correlations between interactions involving residue i
		row = intEnMat[i,:]
		row_mrow = row - row.mean(1)[:,None]
		ssrow = (row_mrow**2).sum(1);
		corrs = np.dot(row_mrow,row_mrow.T)/np.sqrt(np.dot(ssrow[:,None],ssrow[None]))
		sigindices = np.where(corrs > 0.4)
		for m in range(0,len(sigindices[0])):
			row = sigindices[0][m]
			col = sigindices[1][m]
			if row != col and i != row and i != col:
				key = '-'.join(list(map(str,natsorted([i,row])+natsorted([i,col]))))
				if key not in list(sigcorrs.keys()):
					sigcorrs[key] = corrs[row,col]
		for j in range(0,numResidues):
			# Now get all correlations involving residue i and j
			row2 = intEnMat[j,:]
			row2_mrow = row2 - row2.mean(1)[:,None]
			ssrow2 = (row2_mrow**2).sum(1);
			corrs = np.dot(row_mrow,row2_mrow.T)/np.sqrt(np.dot(ssrow[:,None],ssrow2[None]))
			sigindices = np.where(np.abs(corrs) > 0.4)
			for m in range(0,len(sigindices[0])):
				row = sigindices[0][m]
				col = sigindices[1][m]
				if row != col and i != row and j != col and not (i == col and j == row): # Excluding correlations with self.
					key = '-'.join(list(map(str,natsorted([i,row])+natsorted([j,col]))))
					if key not in list(sigcorrs.keys()):
						sigcorrs[key] = corrs[row,col]

		percentCalculated = 25+((i+1)/float(numResidues))*25 # /2 because this is only halfway of calculation.
		logger.info('Interaction energy correlation calculated percentage: %f' % percentCalculated)

		progbar.update()

	logger.info('Calculating interaction energy correlations... completed.')
	logger.info('Collecting significant correlations...')
	df_corr = pandas.DataFrame(columns=['res11','res12','res21','res22','corr'])
	for i in range(0,len(sigcorrs)):
		key = list(sigcorrs.keys())[i]
		matches = re.search('(\d+)-(\d+)-(\d+)-(\d+)',key)
		res11 = int(matches.groups()[0])
		res12 = int(matches.groups()[1])
		res21 = int(matches.groups()[2])
		res22 = int(matches.groups()[3])

		res11_string = getChainResnameResnum(system,res11)
		res12_string = getChainResnameResnum(system,res12)
		res21_string = getChainResnameResnum(system,res21)
		res22_string = getChainResnameResnum(system,res22)

		# Do not include correlations with self.
		#if res11_string == res22_string and res12_string == res21_string:
		#	continue

		df_corr.loc[i] = [res11_string,res12_string,res21_string,res22_string,sigcorrs[key]]
		percentCalculated = 50+(((i+1)/float(len(sigcorrs)))*25)
		logger.info('Interaction energy correlation calculated percentage: %f' % percentCalculated)

	if len(df_corr) == 0:
		logger.info('No significant correlations (with Pearsons r above 0.4) was found. Skipping construction of residue correlation matrix.')
	else:
		logger.info('Saving correlations to file...')
		df_corr.to_csv(outPrefix+'_resIntCorr.csv')
		logger.info('Saving correlations to file... Done.')

		# Constructing the residue correlation matrix
		logger.info('Constructing the residue correlation matrix...')
		rc = np.zeros((numResidues,numResidues))

		for i in range(0,len(sigcorrs)):
			key = list(sigcorrs.keys())[i]
			matches = re.search('(\d+)-(\d+)-(\d+)-(\d+)',key)
			rc_key = np.zeros((numResidues,numResidues))
			if matches:
				res11 = int(matches.groups()[0])
				res12 = int(matches.groups()[1])
				res21 = int(matches.groups()[2])
				res22 = int(matches.groups()[3])

				corr = np.abs(sigcorrs[key])

				rc_key[res11,res21] = corr
				rc_key[res11,res22] = corr
				rc_key[res12,res21] = corr
				rc_key[res12,res22] = corr
				rc_key[res21,res11] = corr
				rc_key[res21,res12] = corr
				rc_key[res22,res11] = corr
				rc_key[res22,res12] = corr

				rc = rc + rc_key
			
			percentCalculated = 75+(((i+1)/float(len(sigcorrs)))*25)
			logger.info('Interaction energy correlation calculated percentage: %f' % percentCalculated)

		logger.info('Constructing the residue correlation matrix... completed.')
		logger.info('Done.')
		np.savetxt(outPrefix+'_resCorr.dat',rc)


def convert_arg_line_to_args(arg_line):
	# To override the same method of the ArgumentParser (to read options from a file)
	# Credit and source: hpaulj from StackOverflow
	# http://stackoverflow.com/questions/29111801/using-fromfile-prefix-chars-with-multiple-arguments-nargs#
	for arg in arg_line.split():
		if not arg.strip():
			continue
		yield arg

if __name__ == '__main__':
	print('Please do not call this method directly. Use grinn -corr instead.')
