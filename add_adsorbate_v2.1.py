# A script to add -O, -OH and -H into Mg surface adsorption sites or subsurface interstitials
# Author: bl4518@ic.ac.uk
# Date: 02/06/2021

import numpy as np
import ase
from ase import Atoms
from ase.build import bulk, hcp0001, add_adsorbate
from ase.io import write
import argparse

parser = argparse.ArgumentParser(description="Building Mg surface with adsorbates")
parser.add_argument('-OH', action='store', dest='OH_numbers', default=0, help='set the number of hydroxyls')
parser.add_argument('-O', action='store', dest='O_numbers', default=0, help='set the number of oxygens')
parser.add_argument('-H', action='store', dest='H_numbers', default=0, help='set the number of hydrogens')
parser.add_argument('-i', action='store', dest='interstitials', help='set the interstitial type')
parser.add_argument('-t', action='store', dest='type', help='cp2k or qe')
args = parser.parse_args()

###################################################################################################

calculation_type = args.type
OH_numbers = int(args.OH_numbers)
O_numbers = int(args.O_numbers)
H_numbers = int(args.H_numbers)
interstitial = args.interstitials

vacuum_length = 10
lattice = dict([('qe_a','3.2094'),('qe_c','5.1862'),('cp2k_a','3.2160'),('cp2k_c','5.2060')])

###################################################################################################
# define a function to build Mg (0001) slab
def slab_builder(layer_numbers, latt_a, latt_c, vacuum_length):
    Mg_0001_slab = hcp0001('Mg',size=(2,2,layer_numbers), a=latt_a, c=latt_c, vacuum=vacuum_length)
    return(Mg_0001_slab)

# define a function to add surface oxygen/hydroxyl
def add_surface_Oxygen(slab, distance, position, offset):
    add_adsorbate(slab, 'O', distance, position, offset)

# define a function to add surface hydroxyls
def add_surface_OH(slab, distance, position, offset):
    bond_length = 0.964
    new_dis = distance + bond_length
    add_adsorbate(slab, 'H', new_dis, position, offset)

# define a function to add H atoms in the tetrahedral-I interstitials
# the distance between the corner atom to the centre in a tetrahedral is about c/2*(0.5+2/3*(a/c)**2)
def add_tetra1_H(slab, latt_a, latt_c, offset):
    bias = - latt_c / 2 * (0.5 - 2 / 3 * (latt_a/latt_c)**2)
    add_adsorbate(slab, 'H', bias, 'hcp', offset)

# define a function to add H atoms in the tetraherdral-II interstitials
def add_tetra2_H(slab, latt_a, latt_c, offset):
    bias = - latt_c / 2 * (0.5 + 2 / 3 * (latt_a/latt_c)**2)
    add_adsorbate(slab, 'H', bias, 'ontop', offset)

# define a function to add H atoms in the octahedral interstitials
def add_octa_H(slab, distance, offset, latt_c):
    new_dis = distance - latt_c / 4
    add_adsorbate(slab, 'H', new_dis, 'fcc', offset)

#################################################################################################
# declare the calculation type: cp2k or quantum-espresso
if calculation_type == 'qe':
    latt_c = float(lattice['qe_c'])
    latt_a = float(lattice['qe_a'])
    slab = slab_builder(8, float(lattice['qe_a']), float(lattice['qe_c']), vacuum_length)
elif calculation_type == 'cp2k':
    latt_c = float(lattice['cp2k_c'])
    latt_a = float(lattice['cp2k_a'])
    slab = slab_builder(8, float(lattice['cp2k_a']), float(lattice['cp2k_c']), vacuum_length)
#################################################################################################
# insert -O/-OH on the surface fcc Hollow sites
offset = ((0,0),(0,1),(1,0),(1,1))

O_total_number = OH_numbers + O_numbers

for i in range(O_total_number):
    add_surface_Oxygen(slab, 0, 'fcc', offset[i])

for i in range(OH_numbers):
    add_surface_OH(slab, 0, 'fcc', offset[i])
################################################################################################
# insert H in the subsurface octahedral/tetrahedral-I/tetrahedral-II interstitials
if interstitial == 'octa':
    for i in range(H_numbers):
        add_octa_H(slab, 0, offset[i], latt_c)

if interstitial == 'tetra1':
    for i in range(H_numbers):
        add_tetra1_H(slab, latt_a, latt_c, offset[i])

if interstitial == 'tetra2':
    for i in range(H_numbers):
        add_tetra2_H(slab, latt_a, latt_c, offset[i])
################################################################################################
# output the structure

write('Mg_0001_OH_H.vasp', slab)
write('Mg_0001_OH_H.xyz', slab)

with open("Mg_0001_OH_H.vasp", "r") as f:
    contents = f.readlines()

value = contents[0]
contents.insert(5, value)

with open("Mg_0001_OH_H.vasp", "w") as f:
    contents = "".join(contents)
    f.write(contents)
