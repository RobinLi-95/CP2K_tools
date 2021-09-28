import os
import sys
import numpy as np
import gen_basis_helpers.analyse_md.traj_core as trajCoreHelp
import gen_basis_helpers.analyse_md.traj_io as trajIoHelp

InpPath = "out_traj.traj"
CurrTraj = trajCoreHelp.readTrajObjFromFileToTrajectoryInMemory(InpPath)
# Dump MD trajactory into Extended-xyz file
trajIoHelp.writeTrajToSimpleExtendedXyzFormat(CurrTraj, "out_traj.exyz")
OutPath = "out_traj.exyz"

file = open(OutPath,'r')
lines = file.readlines()
num_lines = len(lines)

time_series = []
time_line = []

for i in range(num_lines):
    if 'time' in lines[i]:
        time = float(lines[i].split()[-1].split("=")[-1])
        if time > 4500:
            time_series.append(time)
            time_line.append(i)

def StructureWriter(iteration_number):
    n = iteration_number
    element = []
    coord_x = []
    coord_y = []
    coord_z = []
    for i in range(time_line[n]+1,time_line[n]+648):
        element.append(lines[i].split()[0])
        coord_x.append(float(lines[i].split()[1]))
        coord_y.append(float(lines[i].split()[2]))
        coord_z.append(float(lines[i].split()[3]))

    old_number = int(len(element))
################################################################################
# Adding one layer ghost atoms
    fl_z = []
    fl_index = []
    sl_index = []
    sl_x = []
    sl_y = []
    sl_z = []

    aa = np.argsort(coord_z)

    for j in range(72):
        if j < 36:
            fl_z.append(coord_z[aa[j]])
        else:
            sl_x.append(coord_x[aa[j]])
            sl_y.append(coord_y[aa[j]])
            sl_z.append(coord_z[aa[j]])
    tot_z1 = 0
    tot_z2 = 0
    for jj in range(len(fl_z)):
        tot_z1 += fl_z[jj]

    for ii in range(len(sl_z)):
        tot_z2 += sl_z[ii]

    fl_av_z = tot_z1 / len(fl_z)
    sl_av_z = tot_z2 / len(sl_z)
    moving_z = (sl_av_z - fl_av_z) * 2

    for i in range(len(sl_z)):
        element.append('Mg_ghost')
        coord_x.append(sl_x[i])
        coord_y.append(sl_y[i])
        coord_z.append(sl_z[i] - moving_z)

    number_atoms = int(len(element))
################################################################################
    file_name = 'Mg_water_'+str(int(time_series[n]))
    os.system('mkdir %s' % file_name)
    os.system('cp BASIS_MOLOPT %s' % file_name)
    os.system('cp GTH_POTENTIALS %s' % file_name)
    os.system('cp cp2k_batch %s' % file_name)
    os.system('cp template.inp %s' % file_name)

    os.chdir(file_name)
    new_name = file_name + '.inp'
    os.system('mv template.inp %s' % new_name)

    with open('cp2k_batch','r') as f:
        txt = f.readlines()
    out_name = file_name + '.out'
    command = 'mpiexec /rds/general/user/bl4518/home/cp2k-8.2.re/cp2k-8.2/exe/local/cp2k.psmp -i' + ' ' + new_name + ' ' + '-o' + ' ' + out_name
    qsub_name = '#PBS -N' + ' ' + file_name
    txt.insert(2,qsub_name)
    txt.insert(21,command)

    with open('cp2k_batch','w') as f:
        f.writelines(txt)

    output = open('structure.xyz',mode='w')
    output.write("%d\n" % number_atoms)
    output.write("\n")
    for i in range(number_atoms):
        output.write("%s %15.8g %15.8g %15.8g\n" % (element[i],coord_x[i],coord_y[i],coord_z[i]))
    os.chdir('..')

structure_number = len(time_line)

for i in range(structure_number):
    StructureWriter(i)
