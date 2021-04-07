# A script to convert CP2K xyz file into CASTEP cell file with unit cell information
# The print level in CP2K input file should be medium
# ASE is used to read xyz file

import numpy as np
import argparse
import ase
from ase.io import read, write

parser = argparse.ArgumentParser(description='Convert xyz file into cell file')
parser.add_argument('-o', action='store', dest='output_file', help='read cp2k output')
parser.add_argument('-f', action='store', dest='xyz_file', help='read xyz file')
args = parser.parse_args()

outputfile = args.output_file
xyzfile = args.xyz_file

# read the unit cell vectors from CP2K output file
x_vector = []
y_vector = []
z_vector = []
output = open(outputfile,'r')
lines = output.readlines()
num_lines = len(lines)
for i in range(num_lines):
    if 'GENERATE|  Achieved consistency in connectivity generation' in lines[i]:
        x_vector.append(float(lines[i+3].split()[4]))
        x_vector.append(float(lines[i+3].split()[5]))
        x_vector.append(float(lines[i+3].split()[6]))
        y_vector.append(float(lines[i+4].split()[4]))
        y_vector.append(float(lines[i+4].split()[5]))
        y_vector.append(float(lines[i+4].split()[6]))
        z_vector.append(float(lines[i+5].split()[4]))
        z_vector.append(float(lines[i+5].split()[5]))
        z_vector.append(float(lines[i+5].split()[6]))

# read the xyz file with the aid of ASE
xyz = ase.io.read(xyzfile,format='xyz')
positions = xyz.get_positions()
chemical_symbols = xyz.get_chemical_symbols()
atom_number = len(chemical_symbols)

# print out the CASTEP cell file
output = open('structure.cell',mode='w')
output.write("%BLOCK LATTICE_CART\n")
output.write("%15.8g %15.8g %15.8g\n" % (x_vector[0],x_vector[1],x_vector[2]))
output.write("%15.8g %15.8g %15.8g\n" % (y_vector[0],y_vector[1],y_vector[2]))
output.write("%15.8g %15.8g %15.8g\n" % (z_vector[0],z_vector[1],z_vector[2]))
output.write("%ENDBLOCK LATTICE_CART\n")
output.write("\n")
output.write("%BLOCK POSITIONS_ABS\n")
for i in range(atom_number):
    output.write("%s %15.8g %15.8g %15.8g\n" % (chemical_symbols[i],positions[i][0],positions[i][1],positions[i][2]))
output.write("%ENDBLOCK POSITIONS_ABS")
