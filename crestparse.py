#!/usr/bin/python3

"""
crestparse singlefile version v 0.2.2

A command line program to quickly analyze results from crest conformation searches.

Juha Siitonen 14.4.2024
Rice University, Aalto University

TODO:
- Implement multiple simultaneous measurements (currently supports only one at a time)
- Clean up the output formatting functions
- 

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

    # Distance between two cartesian cooardinates in 3D
    def distance(self, atomIndex1, atomIndex2):
        return np.linalg.norm(self.atoms[atomIndex1] - self.atoms[atomIndex2])

    # Angle between three cartesian coordinates in 3D
    def angle(self, atomIndex1, atomIndex2, atomIndex3):
        vector12 = self.atoms[atomIndex1] - self.atoms[atomIndex2]
        vector23 = self.atoms[atomIndex3] - self.atoms[atomIndex2]
        cosineAngle = np.dot(vector12, vector23) / (np.linalg.norm(vector12) * np.linalg.norm(vector23))
        angle = np.arccos(cosineAngle)
        return np.degrees(angle)
    
    # Dihedral angle between four cartesian coordinates in 3D
    def dihedral(self, atomIndex1, atomIndex2, atomIndex3, atomIndex4):
        vector12 = -1*(self.atoms[atomIndex2] - self.atoms[atomIndex1])
        vector32 = self.atoms[atomIndex3] - self.atoms[atomIndex2]
        vector24 = self.atoms[atomIndex4] - self.atoms[atomIndex2]
        vector32 /= np.linalg.norm(vector32)
        projection12 = vector12 - np.dot(vector12, vector32)*vector32
        projection24 = vector24 - np.dot(vector24, vector32)*vector32
        x = np.dot(projection12, projection24)
        y = np.dot(np.cross(vector32, projection12), projection24)
        angle = np.arctan2(y, x)
        return np.degrees(angle)
    
    # Degree of pyramidalization between four cartesian coordinates in 3D
    def pyramidalization(self, atomIndex1, atomIndex2, atomIndex3, atomIndex4):
        angle213 = self.angle(atomIndex2, atomIndex1, atomIndex3)
        angle314 = self.angle(atomIndex3, atomIndex1, atomIndex4)
        angle412 = self.angle(atomIndex4, atomIndex1, atomIndex2)
        totalAngle = angle213 + angle314 + angle412
        # Map the values between a flat plane (360 deg) or tetrahedron (109,5*3 deg)
        return np.interp(totalAngle, [328.5,360], [1,0])

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
    parser.add_argument("-r", "--dihedral", help = "Provides the dihedral in degrees between four atoms with given indeces", nargs=4, type=int)
    parser.add_argument("-p", "--pyramidalization", help = "Provides the pyramidalization degree between four atoms with given indeces. First atom is the central one.", nargs=4, type=int)
    args = parser.parse_args(arguments)

    verbose = args.verbose
    xyzFileIn = args.infile.name
    cutoff = args.cutoff
    temperature = args.temperature
    silent = args.silent
    distances = args.distance
    angles = args.angle
    dihedrals = args.dihedral
    pyramidalization = args.pyramidalization

    conformers = readMultixyzFile(xyzFileIn)
    confomerTotal = len(conformers)

    if verbose:
        logging.basicConfig(level=logging.DEBUG)

    logging.info(f"Successfully read in {confomerTotal} structures")

    calculateRelativeEnergies(conformers)

    logging.info(f"Calculating the Boltzmann distribution at {temperature} K for total of {confomerTotal} conformers...")
    distribution = boltzmannDistribution(conformers, temperature)

    # Apply a cutoff and remove unnesecary conformers
    if cutoff:
        conformers = applyCutoff(conformers, tohartree(float(cutoff)))
        logging.info(f"Applying an energy cutoff of {cutoff} kcal/mol: {len(conformers)} conformers remaining, {confomerTotal-len(conformers)} structures removed")

    # No exports requested?
    if args.extract == None:
        extractionList = []
        
    # Partial list of conformers requested to be exported?
    if args.extract:
        extractionList = args.extract
    
    # All conformers requested to be exported?
    if args.extract == []:
        extractionList = list(range(1, len(conformers)+1))

    # Rough'n'ready implementation, needs some work
    if not silent:
        # Measurement parameter titles
        title = ""
        
        if distances:
            title = f"Distance {distances[0]}-{distances[1]}"
        if angles:
            title = f"Angle {angles[0]}-{angles[1]}-{angles[2]}"
        if dihedrals:
            title = f"Dihedral {dihedrals[0]}-{dihedrals[1]}-{dihedrals[2]}-{dihedrals[3]}"
        if pyramidalization:
            title = f"Pyramidalization {pyramidalization[0]}-{pyramidalization[1]}-{pyramidalization[2]}-{pyramidalization[3]}"
            

        print("#\tE (Hartree)\tdE (Hartree)\tdE (kcal/mol)\tBoltzmann T={0}\t{1}".format(temperature, title))
        
        for i, c in enumerate(conformers):
            # Measurement data
            data = None
            if distances:
                data = c.distance(distances[0], distances[1])
            if angles:
                data = c.angle(angles[0], angles[1], angles[2])
            if dihedrals:
                data = c.dihedral(dihedrals[0], dihedrals[1], dihedrals[2], dihedrals[3])
            if pyramidalization:
                data = c.pyramidalization(pyramidalization[0], pyramidalization[1], pyramidalization[2], pyramidalization[3])
            
            if data:
                print('{0:d}\t{1:f}\t{2:f}\t{3:f}\t{4:%}\t\t{5:f}'.format(
                    c.index,
                    c.energy,
                    c.relativeEnergy,
                    tokcal(c.relativeEnergy),
                    distribution[i],
                    data))
            else:
                print('{0:d}\t{1:f}\t{2:f}\t{3:f}\t{4:%}'.format(
                    c.index,
                    c.energy,
                    c.relativeEnergy,
                    tokcal(c.relativeEnergy),
                    distribution[i]))

            
    if extractionList:
        logging.info("Exporting files...")
        
        for i in extractionList:
            path = f"conf_{i}.xyz"
            writexyzFile(conformers[i-1], path)
            if not silent:
                print("Conformer #{0} exported to file {1}".format(i, path))

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
