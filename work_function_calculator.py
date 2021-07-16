import numpy as np
from scipy import interpolate
import sys
import os
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='A script to calculate the work function in direct/macroscopic way')
parser.add_argument('-f', action='store', dest='filename', help='declare the filename')
parser.add_argument('-t', action='store', dest='type', help='declare the type of functionals')
parser.add_argument('-e', action='store', dest='fermi_energy', help='declare the fermi energy')
args = parser.parse_args()

input = args.filename
functional = args.type
###################################################################################################################
# interolating more points in the original electrostatic potential data

data = np.loadtxt(input)
x = data[:,0]
y = data[:,1]

length = len(x)
xnew = np.linspace(0,x[length-1],num=5000)

f = interpolate.interp1d(x,y,kind='cubic')
ynew = f(xnew)

scale = len(xnew)

outputfile = open("interpolate.txt", mode='w')
outputfile.write("#  Distance(Ang)     Plannar_Potential(ha)\n")
for i in range(0,scale,1):
    outputfile.write("%15.8g %15.8g\n" % (xnew[i],ynew[i]))
outputfile.close()

###################################################################################################################
# calculating the macroscopic average electrostatic potential

# define some important value
data = numpy.loadtxt('interpolate.txt')
distance = data[:,0]
pot = data[:,1] * 27.2113961318
grid_interval = data[1,0]
number_of_points = int(average_window/grid_interval/2)

if functional == 'PBE':
    average_window = 5.206
elif functional == 'BLYP':
    average_window = 5.124

length = len(distance) - 1
finaldistance = distance[length]
TopOfSlab = finaldistance - 10
EndOfSlab = 10

#find the integral range
for i in range(0,length,1):
    if distance[i] < int(TopOfSlab + 1) and distance[i] > int(TopOfSlab - 1):
        end_index = i

for i in range(0,length,1):
    if distance[i] < 10.1 and distance[i] > 10:
        initial_index = i

half_slab_length = (distance[end_index] - distance[initial_index]) * 0.5
inital_value = distance[initial_index] + half_slab_length * 0.5
end_value = distance[end_index] - half_slab_length * 0.5

for ii in range(0,length,1):
    if distance[ii] < int(inital_value + 1) and distance[ii] > int(inital_value - 1):
        integer_bottom = ii

for iii in range(0,length,1):
    if distance[iii] < int(end_value + 1) and distance[iii] > int(end_value - 1):
        integer_top = iii
inner_distance = integer_top - integer_bottom + 1

sum_pot = 0
macro_pot = numpy.zeros(inner_distance,dtype=float)

#do the macrocopic average
for a in range(integer_bottom,integer_top+1,1):
    new_bottom_index = a - number_of_points
    new_top_index = a + number_of_points
    sum_pot = 0
    for b in range(new_bottom_index,new_top_index+1,1):
        sum_pot = sum_pot + grid_interval * pot[b]
    macro_pot[a-integer_bottom] = sum_pot / average_window

#output the file
z_direction = numpy.zeros(inner_distance,dtype=float)
for i in range(0,inner_distance,1):
    z_direction[i] = distance[integer_bottom+i]

outputfile = open("macroscopic_new.txt", mode='w')
outputfile.write("#  Distance(Ang)     Macrocopic_Potential(ha)\n")
for i in range(0,inner_distance,1):
    outputfile.write("%15.8g %15.8g\n" % (z_direction[i],macro_pot[i]))
outputfile.close()

###################################################################################################################
# plot the electrostatic potential
data_2 = np.loadtxt('macroscopic_new.txt')


ha = 27.2113961318
distance = data[:,0]
pot = data[:,1] * ha
fermi_slab = float(args.fermi_energy) * ha

plt.title('Electrostatic potential of Mg (0001) surface')
plt.plot(distance, pot, color='blue', label='Average Potential in XY Plane')
plt.plot(data_2[:,0],data_2[:,1],color='orange',label='Macroscopic Average Potential')
plt.xlabel('Distance_Z(Angstrom)')
plt.ylabel('Electrostatic Potential(eV)')
plt.hlines(Fermi_energy, 0, distance[len(distance) - 1], colors="r",linestyles ="dashed",label='Fermi Energy')
plt.xlim(0,distance[len(distance) - 1])
plt.legend()
plt.savefig('electrostatic_mg_basal',dpi=500)
plt.show()

###################################################################################################################
# calculate the work function
# calculate the macroscopic potential inside the slab
length = len(data_2[:,0])
tol = 0

for i in range(0,length,1):
    tol = tol + data_2[i,1]

macro_average_slab = tol / length

# calculate the vacuum potential
pot_vacuum = pot[0]

if functional == 'PBE':
    macro_slab_bulk = 0.60067999
    fermi_bulk = 0.10632975843621 * ha
elif functional == 'BLYP':
    macro_slab_bulk = 0.58271422
    fermi_bulk = 0.12032507110167 * ha

# calculate the work function in two methods
direct_work_function = pot_vacuum - fermi_slab
corrected_work_function = pot_vacuum - macro_average_slab + macro_average_bulk - fermi_bulk

print('The direct work function is: %15.8g' % direct_work_function )
print('The correct work function is: %15.8g' % corrected_work_function )
