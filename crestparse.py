#!/usr/bin/python3

"""
crestparse singlefile version v 0.2

A command line program to quickly analyze results from crest conformation searches.

Juha Siitonen 13.4.2024
Rice University, Aalto University

This project is licensed under the terms of the MIT license.
"""

import os
import sys
import argparse
import re
import logging
import math
import numpy as np


class Conformer:
    def __init__(self, index, energy, xyzFile):
        self.index = index
        self.energy = energy
        self.relativeEnergy = 0
        self.xyzFile = xyzFile
        self.boltzmannFactor = 0
        self.atoms = []
        
        # Populate the atom table, skip the two first rows
        for row in xyzFile[2:]:
            atomData = row.split()
            
            x = float(atomData[1])
            y = float(atomData[2])
            z = float(atomData[3])
			
            self.atoms.append(np.array([x,y,z]))
        
    def formatxyz(self):
        stringOutput = ""
        for line in self.xyzFile:
            stringOutput += line + "\n"
        return stringOutput

    def distance(self, atomIndex1, atomIndex2):
        return np.linalg.norm(self.atoms[atomIndex1] - self.atoms[atomIndex2])

    def angle(self, atomIndex1, atomIndex2, atomIndex3):
        vector12 = self.atoms[atomIndex1] - self.atoms[atomIndex2]
        vector23 = self.atoms[atomIndex3] - self.atoms[atomIndex2]
        cosineAngle = np.dot(vector12, vector23) / (np.linalg.norm(vector12) * np.linalg.norm(vector23))
        angle = np.arccos(cosineAngle)
        return np.degrees(angle)
    

# Auxiliary functions for energies
def tokcal(e):
    return e*627.5

def tohartree(e):
    return e/627.5

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
    parser.add_argument("-e", "--extract", help = "Conformer indexes to be extracted into separate files (conf_n.xyz). If not provided, all conformers will be extracted.", nargs="*", type=int)
    parser.add_argument("-d", "--distance", help = "Provides the distance in angstrom between two atoms with given indeces", nargs=2, type=int)
    parser.add_argument("-a", "--angle", help = "Provides the angle in degrees between three atoms with given indeces", nargs=3, type=int)
    args = parser.parse_args(arguments)

    verbose = args.verbose
    xyzFileIn = args.infile.name
    cutoff = args.cutoff
    temperature = args.temperature
    silent = args.silent
    distances = args.distance
    angles = args.angle

    conformers = readMultixyzFile(xyzFileIn)
    confomerTotal = len(conformers)

    if verbose:
        logging.basicConfig(level=logging.DEBUG)

    logging.info("Successfully read in " + str(confomerTotal) + " structures")

    calculateRelativeEnergies(conformers)

    logging.info("Calculating the Boltzmann distribution at " + str(temperature) + " K for total of "+ str(confomerTotal) + " conformers...")
    distribution = boltzmannDistribution(conformers, temperature)

    # Apply a cutoff and remove unnesecary conformers
    if cutoff:
        conformers = applyCutoff(conformers, tohartree(float(cutoff)))
        logging.info("Applying an energy cutoff of " + str(cutoff) + " kcal/mol: " + str(len(conformers)) + " conformers remaining, " + str(confomerTotal-len(conformers)) + " structures removed")

    # No extraction requested?
    if args.extract == None:
        extractionList = []
        
    # Partial list requested?
    if args.extract:
        extractionList = args.extract
    
    # All requested?
    if args.extract == []:
        extractionList = list(range(1, len(conformers)+1))

    # Rough'n'ready implementation, needs some work
    if not silent:
        if distances:
            print("#\tE (Hartree)\tdE (Hartree)\tdE (kcal/mol)\tBoltzmann T=" + str(temperature) + " K\tDistance " + str(distances[0]) + "-" + str(distances[1]))
        elif angles:
            print("#\tE (Hartree)\tdE (Hartree)\tdE (kcal/mol)\tBoltzmann T=" + str(temperature) + " K\tAngle " + str(angles[0]) + "-" + str(angles[1]) + "-" + str(angles[2]))
        else:
            print("#\tE (Hartree)\tdE (Hartree)\tdE (kcal/mol)\tBoltzmann T=" + str(temperature) + " K")
        for i, c in enumerate(conformers):
            if distances:
                print('{0:2d}\t{1:f}\t{2:f}\t{3:f}\t{4:%}\t\t{5:f}'.format(c.index, c.energy, c.relativeEnergy, tokcal(c.relativeEnergy), distribution[i], c.distance(distances[0], distances[1])))
            elif angles:
                print('{0:2d}\t{1:f}\t{2:f}\t{3:f}\t{4:%}\t\t{5:f}'.format(c.index, c.energy, c.relativeEnergy, tokcal(c.relativeEnergy), distribution[i], c.angle(angles[0], angles[1], angles[2])))
            else:
                print('{0:2d}\t{1:f}\t{2:f}\t{3:f}\t{4:%}'.format(c.index, c.energy, c.relativeEnergy, tokcal(c.relativeEnergy), distribution[i]))

    logging.info("Exporting files...")
    if extractionList:
        for i in extractionList:
            path = "conf_" + str(i) + ".xyz"
            writexyzFile(conformers[i-1], path)
            if not silent:
                print("Conformer #" + str(i) + " exported to file " + path)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
