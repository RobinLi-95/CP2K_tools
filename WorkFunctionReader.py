import numpy as np
import os
import sys
import re
import subprocess
import matplotlib.pyplot as plt
from ase.units import Bohr, Ry, Hartree
from ase.io import read, write
from ase.io.cube import read_cube
from scipy.stats import norm
import time
import datetime

def CubeFileCheck():
    CheckIndex = 0
    FileList = os.listdir()

    if 'Mg_0001-v_hartree-1_0.cube' in FileList:
        CheckIndex = 1
    else:
        print('CubeFile is not in the current directory.')

    return(CheckIndex)

def CubeFileReader(CubeFile, index):
    file = open(CubeFile)
    tot_potential = read_cube(file,read_data=True)
    potl = tot_potential['data']
    atoms = tot_potential['atoms']

    # Read in lattice parameters and scale factor
    cell = atoms.cell

    lattice = np.dot(cell, cell.T).diagonal()
    lattice_length = lattice ** 0.5

    # Read in potential data
    ngridpts = np.array(potl.shape)
    totgridpts = ngridpts.prod()

    # Perform average
    a = 0
    b = 1
    idir = 2

    a = (idir + 1) % 3
    b = (idir + 2) % 3

    average = np.zeros(ngridpts[idir])
    for ipt in range(ngridpts[idir]):
        average[ipt] = potl[:,:,ipt].sum()

    average /= ngridpts[a] * ngridpts[b]
    converted_avg = average * Hartree

    # Print out the average
    filename = 'Potential_' + str(index) + '.txt'
    outputfile = open(filename,"w")
    outputfile.write("# Distance(Ang)   Potential(eV)\n")
    xdiff = lattice_length[idir]/float(ngridpts[idir])

    for i in range(ngridpts[idir]):
        x = i * xdiff
        outputfile.write("%15.8g %15.8g\n" % (x,converted_avg[i]))
    outputfile.close()

    return()

def PotentialPlot(filename, Fermi, index):
    data = np.loadtxt(filename)

    z_dir = data[:,0]
    pot = data[:,1]

    z_lim = z_dir[len(z_dir) - 1]
    title = 'Potential Plot of Mg/Water at' + str(index) + ' ' +'ps'
    header = 'Potential_at_' + str(index)

    plt.title(title)
    plt.plot(z_dir, pot, color='blue', label='Average Potential in XY Plane')
    plt.hlines(Fermi, 0, z_lim, colors='r', linestyle="dashed", label='Fermi Energy')
    plt.xlim(0,z_lim)

    plt.xlabel('Z(Angstrom)')
    plt.ylabel('Potential(eV)')
    plt.legend()
    plt.savefig(header, dpi=500)
    plt.close()

    return()

def FermiCalc(CP2KOutFile):
    fermi = 0
    with open(CP2KOutFile,'r') as file:
        for line in file:
            if re.search('Fermi', line):
                fermi = float(line.split()[2]) * Hartree

    return(fermi)

def OutputFileCheck(CurrentDir):
    CheckIndex = 0
    FileList = os.listdir()
    OutPut = CurrentDir + '.out'

    if OutPut in FileList:
        CheckIndex = 1
    else:
        print('%s is not in the current directory.' % OutPut)

    return(CheckIndex)

def WorkFuctionCalc(filename, fermi_level):
    data = np.loadtxt(filename)

    z_dir = data[:,0]
    pot = data[:,1]

    tot_pot = 0
    num = 0

    for i in range(len(z_dir)):
        if z_dir[i] > 42 and z_dir[i] < 48:
            tot_pot += pot[i]
            num += 1

    vac_pot = tot_pot / num
    PZC = vac_pot - fermi_level

    return(PZC)

def FileWriter(xcol, ycol):
    filename = 'PZC_MD.txt'
    mean_pzc = np.mean(ycol)
    output = open(filename,mode='w')
    output.write('# Potential at Zero Charge of Mg/Water\n')
    output.write('The average potential of zero charge is %15.8g\n' % mean_pzc)
    output.write('\n')

    for i in range(len(xcol)):
        output.write('%s %15.8g\n' % (str(xcol[i]),ycol[i]))

    return()

def GeneralPlot(timestep, pzc):
    plt.title('Potential of Zero Charge for Mg/Water')
    plt.plot(timestep, pzc, color='r', label='Potential of Zero Charge for Mg Basal Plane')
    plt.xlabel('Time Step (fs)')
    plt.ylabel('PZC (eV)')

    plt.legend()
    plt.savefig('pzc_vs_time.png', dpi=400)
    plt.close()

    return()

def NormDistributedFit(data):
    std = np.std(data, ddof = 1)
    mean = np.mean(data)

    domain = np.linspace(np.min(data),np.max(data))
    plt.plot(domain, norm.pdf(domain,mean,std),
             label='$\mathcal{N}$' + f'$( \mu \\approx {round(mean)}, \sigma \\approx {round(std)})$')
    plt.hist(data, edgecolor='black', alpha=.5, density=True)
    plt.title('Normal Fit of PZC')
    plt.xlabel('Potential of Zero Charge(eV)')
    plt.ylabel('Density')
    plt.legend()
    plt.savefig('Potential_zero_charge.png',dpi=500)
    plt.close()

    return()


#########################################################################################################
# Main part of this code
starttime = time.perf_counter()
print ("Starting calculation at")
print (time.strftime("%H:%M:%S on %a %d %b %Y"))

Timestep = []
PotentialZeroCharge = []
Template = 'Mg_water_'

for i in range(4536,20536,100):
    CurrDir = Template + str(i)
    os.chdir(CurrDir)
    output = CurrDir + '.out'

    outindex = OutputFileCheck(CurrDir)
    cubeindex = CubeFileCheck()
    fermi = FermiCalc(output)

    if outindex == 0 or cubeindex == 0 or fermi == 0:
        os.chdir('..')
    elif outindex == 1 and cubeindex == 1 and fermi != 0:
        print('Current directory is %s' % CurrDir)
        Timestep.append(i)
        CubeFileReader('Mg_0001-v_hartree-1_0.cube', i)
        PotentialFile = 'Potential_' + str(i) + '.txt'

        PotentialZeroCharge.append(WorkFuctionCalc(PotentialFile,fermi))
        PotentialPlot(PotentialFile,fermi,i)

        os.chdir('..')

FileWriter(Timestep, PotentialZeroCharge)
GeneralPlot(Timestep, PotentialZeroCharge)
NormDistributedFit(PotentialZeroCharge)

print ("Calculation is done.")

endtime = time.perf_counter()
runtime = endtime - starttime

print ("\nEnd of calculation.")
print ("Program was running for %.2f seconds." % runtime)
