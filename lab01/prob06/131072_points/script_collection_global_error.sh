#!/bin/bash

index_min=1
index_max=17
delta_index=1
pot=2

for i in $(seq ${index_min} ${delta_index} ${index_max});
	do
		step=$(echo "${pot}^${i}" | bc)
		array='-n '${step},1p''
		sed ${array} result_y1.dat >> result_global_error.dat
	done

# Notes:

# compilation:
# 	1) chmod +x script.sh
#	2) ./script.sh
