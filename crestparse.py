#!/usr/bin/python3

"""
crestparse singlefile version v 0.1

A command line program to quickly analyze results from crest conformation searches.

Juha Siitonen 2.5.2020
Rice University

This project is licensed under the terms of the MIT license.
"""

import os
import sys
import argparse
import re
import logging
import math

def tokcal(e):
	return e*627.5

def tohartree(e):
	return e/627.5

class Conformer:
	def __init__(self, index, energy, xyzFile):
		self.index = index
		self.energy = energy
		self.relativeEnergy = 0
		self.xyzFile = xyzFile
		self.boltzmannFactor = 0

	def formatxyz(self):
		stringOutput = ""
		for line in self.xyzFile:
			stringOutput += line + "\n"
		return stringOutput

def getMinimum(confList):
	return min(confList, key=lambda c: c.energy)

def conformerEnergyDifference(conf1, conf2):
	return conf1.energy - conf2.energy

def calculateRelativeEnergies(confList):
	lowestConformer = getMinimum(confList)
	for c in confList:
		c.relativeEnergy = conformerEnergyDifference(c, lowestConformer)

def boltzmannDistribution(confList, temperature):
	R = 0.001987204

	totalEnergy = sum(tokcal(c.relativeEnergy) for c in confList)
	factorList = [math.exp(-tokcal(c.relativeEnergy/(R*temperature))) for c in confList]
	factorSum = sum(factorList)
	distribution = [factor/factorSum for factor in factorList]

	return distribution

def applyCutoff(confList, energyCutoff):
	return [conformer for conformer in confList if conformer.relativeEnergy < energyCutoff]

def readMultixyzFile(filename):
	with open(filename) as f:
		lines = f.read().splitlines()

	conformerList = []
	index = 1
	dataLength = int(lines[0]) + 2

	for structureRow in range(0,len(lines),dataLength):
		conformerList.append(Conformer(int(index), float(lines[structureRow+1]), lines[structureRow:structureRow+dataLength]))
		index += 1

	return conformerList

def writexyzFile(conformer, filename):
	with open(filename, "w") as f:
		f.write(conformer.formatxyz())

class commaSeparateAction(argparse.Action):
	def __call__(self, parser, namespace, values, option_string = None):
		indexes = values.split(",")
		setattr(namespace, self.dest, [int(i) for i in indexes])

def main(arguments):
	parser = argparse.ArgumentParser(description = __doc__, formatter_class = argparse.RawDescriptionHelpFormatter)
	parser.add_argument("infile", help = "Multistructure xyz file from crest (typically crest_conformers.xyz)", type = argparse.FileType('r'))
	parser.add_argument("-v", "--verbose", help = "Verbose mode to display additional details", action = "store_true")
	parser.add_argument("-s", "--silent", help = "Run in silent mode, no table output", action = "store_true")
	parser.add_argument("-c", "--cutoff", help = "Add an energy cutoff in kcal/mol, only conformers up to the cutoff will be shown", type = float)
	parser.add_argument("-t", "--temperature", help = "Set the Boltzmann distribution temperature (in K)", type = float, default = 298.15)
	parser.add_argument("-e", "--extract", help = "Conformer indexes to be extracted into separate files (conf_n.xyz). Either 4 or 4,3,6", type = str, action = commaSeparateAction)
	args = parser.parse_args(arguments)
	
	extractionList = args.extract
	verbose = args.verbose
	xyzFileIn = args.infile.name
	cutoff = args.cutoff
	temperature = args.temperature
	silent = args.silent

	conformers = readMultixyzFile(xyzFileIn)
	confomerTotal = len(conformers)

	if verbose:
		logging.basicConfig(level=logging.DEBUG)
	
	logging.info("Successfully read in " + str(confomerTotal) + " structures")

	calculateRelativeEnergies(conformers)

	logging.info("Calculating the Boltzmann distribution at " + str(temperature) + " K for total of "+ str(confomerTotal) + " conformers...")
	distribution = boltzmannDistribution(conformers, temperature)

	if cutoff:
		conformers = applyCutoff(conformers, tohartree(float(cutoff)))
		logging.info("Applying an energy cutoff of " + str(cutoff) + " kcal/mol: " + str(len(conformers)) + " conformers remaining, " + str(confomerTotal-len(conformers)) + " structures removed")

	if not silent:
		print("#\tE (Hartree)\tdE (Hartree)\tdE (kcal/mol)\tBoltzmann (%) T = " + str(temperature) + " K")
		for i, c in enumerate(conformers):
			print('{0:2d}\t{1:f}\t{2:f}\t{3:f}\t{4:%}'.format(c.index, c.energy, c.relativeEnergy, tokcal(c.relativeEnergy), distribution[i]))

	logging.info("Exporting files...")
	if extractionList:
		for i in extractionList:
			path = "conf_" + str(i) + ".xyz"
			writexyzFile(conformers[i-1], path)
			print("Conformer #" + str(i) + " exported to file " + path)

if __name__ == '__main__':
	sys.exit(main(sys.argv[1:]))
