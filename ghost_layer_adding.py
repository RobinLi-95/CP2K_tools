# This is a script to construct the Mg (0001) slab with different numbers of ghost atom layers

import numpy as np
import argparse
import sys
import os

parser = argparse.ArgumentParser(description='A script to adding the ghost layers' )
parser.add_argument('-f', action="store", dest="filename", help='declare the filename')
args = parser.parse_args()

# read the input file
input = args.filename
f = open(input,"r")
lines = f.readlines()
length_of_file = len(lines)

# determine the coordinates part in the input file
top_layer_x = []
top_layer_y = []
top_layer_z = []

bottom_layer_x = []
bottom_layer_y = []
bottom_layer_z = []
for i in range(0,length_of_file):
    if (len(lines[i]) !=0) and ('COORD' in lines[i]) and ('END' not in lines[i]):
        words_1layer = lines[i+1].split()
        words_2layer = lines[i+2].split()
        bottom_layer_x.append(float(words_1layer[1]))
        bottom_layer_x.append(float(words_2layer[1]))
        bottom_layer_y.append(float(words_1layer[2]))
        bottom_layer_y.append(float(words_2layer[2]))
        bottom_layer_z.append(float(words_1layer[3]))
        bottom_layer_z.append(float(words_2layer[3]))
    if (len(lines[i]) !=0 and ('END COORD' in lines[i])):
        words_3layer = lines[i-1].split()
        words_4layer = lines[i-2].split()
        top_layer_x.append(float(words_3layer[1]))
        top_layer_x.append(float(words_4layer[1]))
        top_layer_y.append(float(words_3layer[2]))
        top_layer_y.append(float(words_4layer[2]))
        top_layer_z.append(float(words_3layer[3]))
        top_layer_z.append(float(words_4layer[3]))

c_axis_distance = bottom_layer_z[1] - bottom_layer_z[0]


# building ghost layers
ghost_kind = 'Mg_ghost'
ghost_x_top = []
ghost_y_top = []
ghost_z_top = []
ghost_x_bottom = []
ghost_y_bottom = []
ghost_z_bottom = []
for n_ghost in range(1,4,1):
    if (n_ghost % 2 != 0):
        ghost_x_top.append(top_layer_x[1])
        ghost_y_top.append(top_layer_y[1])
        ghost_z_top.append(top_layer_z[0] + n_ghost * c_axis_distance)

        ghost_x_bottom.append(bottom_layer_x[1])
        ghost_y_bottom.append(bottom_layer_y[1])
        ghost_z_bottom.append(bottom_layer_z[0] - n_ghost * c_axis_distance)
    else:
        ghost_x_top.append(top_layer_x[0])
        ghost_y_top.append(top_layer_y[0])
        ghost_z_top.append(top_layer_z[0] + n_ghost * c_axis_distance)

        ghost_x_bottom.append(bottom_layer_x[0])
        ghost_y_bottom.append(bottom_layer_y[0])
        ghost_z_bottom.append(bottom_layer_z[0] - n_ghost * c_axis_distance)

# write new input file
for i in range(3):
    os.system('cp Mg_0001.inp Mg_0001_'+str(i+1)+'.inp')

for i in range(3):
    with open('Mg_0001_'+str(i+1)+'.inp','r') as f:
        a = f.readlines()
    for ilayer in range(i+1):
        # print(ilayer)
        atom_layer_top = ghost_kind + ' ' + str(ghost_x_top[ilayer]) + ' '+ str(ghost_y_top[ilayer]) + ' '+ str(ghost_z_top[ilayer]) + '\n'
        a.insert(20, atom_layer_top)
        atom_layer_bottom = ghost_kind + ' ' + str(ghost_x_bottom[ilayer]) + ' '+ str(ghost_y_bottom[ilayer]) + ' '+ str(ghost_z_bottom[ilayer]) + '\n'
        a.insert(20, atom_layer_bottom)
# insert the ghost atom kind part in the inout file
    kind = ' ' + ' ' + ' ' + ' ' + '&KIND Mg_ghost' + '\n'
    element = ' ' + ' ' + ' ' + ' ' + ' ' + 'ELEMENT  Mg' + '\n'
    basis = ' ' + ' ' + ' ' + ' ' + ' ' + 'BASIS_SET rc-7pt5-sspd+s-single-exp-atom-opt' + '\n'
    ghost = ' ' + ' ' + ' ' + ' ' + ' ' + 'GHOST .TRUE.' + '\n'
    end_kind = ' ' + ' ' + ' ' + ' ' + '&END KIND' + '\n'
    a.insert(14,kind)
    a.insert(15,element)
    a.insert(16,basis)
    a.insert(17,ghost)
    a.insert(18,end_kind)

    with open('Mg_0001_'+str(i+1)+'.inp','w') as f:
        f.writelines(a)
