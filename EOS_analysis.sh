#!/bin/bash

cutoff=300
ref_lattice_a=3.21
interval=0.01

input_file=brucite.in
output_file=brucite.out
plot_file=lattice_data.txt


echo "# lattice_a vs total energy" > $plot_file
echo "# Date: $(date)" >> $plot_file
echo "# PWD: $PWD" >> $plot_file
echo "# CUTOFF: ${cutoff}"  >> $plot_file
echo -n "# Lattice_a (Angstrom) | Total Energy (Ry)" >> $plot_file
echo "\n" >> $plot_file

for ii in {-10..10..1} ; do
    aa=$(expr ${ref_lattice_a}+${ii}*${interval} |  bc)
    work_dir=a_${aa}
    total_energy=$(grep -e '!' $work_dir/$output_file | awk '{print $5}')
    printf "%10.2f  %15.10f" $ii $total_energy >> $plot_file
    printf "\n" >> $plot_file
done
