import sys, numpy, math
import argparse
from scipy.optimize import leastsq

##################################
Ha = 27.211396641308
Ry = 13.605698320654

parser = argparse.ArgumentParser(description='A script to fit eos of hexagonal/tetragonal crystal')
parser.add_argument('-t', action='store', dest='l_type', help='declare the type of crystal, could be hexagonal or tetragonal')
parser.add_argument('-r', action='store', dest='c_a_ratio', help='declare the c_over_a ratio')
args = parser.parse_args()

##############################################################################################################################
# Birch-Murnaghan equation of state
def eos_birch_murnaghan(params, vol):
    'From Phys. Rev. B 70, 224107'
    E0, B0, Bp, V0 = params
    eta = (V0/vol)**(1.0/3.0)
    E = E0 + 9.0*B0*V0/16.0 * ((eta**2-1.0)**3 * Bp + (eta**2-1.0)**2 * (6.0 - 4.0*eta**2))
    return E


def hexagonal_vol_calc(a, c_ratio):
    volume = math.sqrt(3) * c_ratio * (a ** 3) / 2 / 2
    return volume

def tetragonal_vol_calc(a, c_ratio):
    volume = c_ratio * (a ** 3)
    return volume
  
def print_params(label, params):
    E0, B0, Bp, V0 = params
    print(label, ": E0 = %f eV" % (E0))
    print(label, ": B0 = %f GPa" % (B0*160.21765))
    print(label, ": Bp = %f" % (Bp))
    print(label, ": V0 = %f angstrom^3" % (V0))
    print

################################################################################################
# transform the lattice parameter to unit cell volume
vol = []
data = numpy.loadtxt('lattice_data.txt')
lattice_a = numpy.array(data[:,0])
ene = numpy.array(data[:,1]) * Ha
lattice_type = args.l_type
ratio = args.c_a_ratio

if lattice_type=='hexagonal':
  for i in range(len(lattice_a)):
    vol.append(hexagonal_vol_calc(lattice_a[i],ratio))
elif lattice_type == 'tetragonal':
  for i in range(len(lattice_a)):
    vol.append(tetragonal_vol_calc(lattice_a[i],ratio))

################################################################################################
    
# fit a parabola to the data and get inital guess for equilibirum volume and bulk modulus
# note the bulk modulus here referes to the linear moduli of hexagonal/tetragonal along a-aixs (C11 in the elastic matrix)
a, b, c = numpy.polyfit(vol, ene, 2)
V0 = -b/(2*a)
E0 = a*V0**2 + b*V0 + c
B0 = 2*a*V0
Bp = 4.0

# initial guesses in the same order used in the Murnaghan function
x0 = [E0, B0, Bp, V0]


# fit the equations of state
target = lambda params, y, x: y - eos_birch_murnaghan(params, x)
birch_murn, ier = leastsq(target, x0, args=(ene,vol))
print_params("Birch-Murnaghan", birch_murn)

#################################################################################################
# plot the fitting curve
import pylab
vfit = numpy.linspace(min(vol),max(vol),100)

pylab.plot(vol, ene, 'ro')
pylab.plot(vfit, eos_birch_murnaghan(birch_murn,vfit), label='Birch-Murnaghan',color = 'blue')
pylab.xlabel('Volume per Atom ($\AA^3$)')
pylab.ylabel('Energy per Atom (eV)')
pylab.legend(loc='best',fontsize='large')
pylab.savefig('eos_bm_fitting',dpi=600)
pylab.show()
