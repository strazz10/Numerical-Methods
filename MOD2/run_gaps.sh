#!/bin/bash	
export LC_NUMERIC="en_US.UTF-8"      #sistema le virgole dei decimali
                                     #i file vanno ottenuti compilando almeno una volta compiler.sh

gfortran gaps.f90 davidson.guess.lib2.o lapack_tebd.lib.f -O2 -o gaps.x   #compila linkando i file necessari

#misure a L fissato
start=0.005
end=0.050
spacing=0.005

output_file="gaps_5.txt"                         #mettere L prima di ogni run (gaps_L.txt), segnarsi se obc, pbc
: > "$output_file"                                #pulisce il file
                          
for g in $(seq $start $spacing $end); do          #loop per i valori del campo g desiderati 
    echo "Sending input g= $g to simulation" 
    echo "$g" | ./gaps.x >> "$output_file"        #pipe
done
echo "Done simulating"

#misure a g fissato 
#start=3
#end=22
#spacing=1

#output_file="gaps_1.25_d0_pbc.dat"                   #segnarsi il valore di g e se obc, pbc
#: > "$output_file"                                #pulisce il file
                          
#for L in $(seq $start $spacing $end); do          #loop per i valori della lunghezza desiderati 
#    echo "Sending input L= $L to simulation" 
#    echo "$L" | ./gaps.x >> "$output_file"        #pipe
#done
#echo "Done simulating"
