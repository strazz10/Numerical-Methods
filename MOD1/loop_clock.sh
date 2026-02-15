#!/bin/bash	
export LC_NUMERIC="en_US.UTF-8"      #sistema le virgole dei decimali

gfortran clock.f90 -O2 -o clock.x

#loop through desired numbers 
for i in $(seq 0.80 0.01 0.94); do
    output_file="clock_beta_${i}_16.txt"     #mettere L prima di ogni run
    echo "Sending input beta= $i to simulation"
    # pipe 
    echo "$i" | ./clock.x > "$output_file"
done
echo "Done simulating"

