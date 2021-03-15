#!/bin/bash
# This is a script for generating a series CP2K input file to obtain the Volume-Energy curves for Hexagonal/Tetragonal crystals

# Define the reference lattice parameter a and c/a ratio for your structure 
reference_a=3.21
interval=0.01
ratio=1.62

# Define the basis set/pseudopotential/template input/batch file for your calculation 
basis_file=BASIS_MOLOPT
potential_file=GTH_POTENTIALS
template_file=template.inp
input_file=Mg.inp
batch_file=cp2k_batch

# Generate the input file by varing the a/c parameters
for ii in {-10..10..1} ; do
    aa=$(expr ${reference_a}+${interval}*${ii} | bc)
    cc=$(expr ${aa}*${ratio} | bc)
    work_dir=a_${aa}
    if [ ! -d $work_dir ] ; then
       mkdir $work_dir
    else 
       rm -r $work_dir/*
    fi
    sed -e "s/AA/${aa}/g" \
    sed -e "s/CC/${cc}/g" \
        $template_file > $work_dir/$input_file
    cp $basis_file $work_dir
    cp $potential_file $work_dir
    cp $batch_file $work_dir
done
