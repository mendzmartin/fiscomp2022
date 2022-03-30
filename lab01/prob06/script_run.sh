#!/bin/bash

# module paths
	pth_mod_01='../../modules/module_presition.f90'
	pth_mod_02='../../modules/module_functions_1D.f90'
	pth_mod_03='../../modules/module_functions_2D.f90'
	pth_mod_04='../../modules/module_EDO_segundo_orden.f90'
	pth_mod_05='../../modules/module_numerical_error.f90'

	pth_mod=${pth_mod_01}' '${pth_mod_02}' '${pth_mod_03}' '${pth_mod_04}' '${pth_mod_05}

# object code name
	ob_cod_name='ODE_second_order_v2.o'

# fortran code name
	f90_cod_name='ODE_second_order_v2.f90'

#flags compiles
#{see: https://faculty.washington.edu/rjl/classes/am583s2013/notes/gfortran_flags.html}
	## warning flags
		flg_w01='-Wall -Wextra -Wconversion -pedantic'
	## debugging flags
		flg_d01='-ffpe-trap=zero -ffpe-trap=overflow -ffpe-trap=underflow' 
		flg_d02='-g -fbacktrace -fbounds-check'
	## optimization flags
		flg_o01='-o3 -ftree-vectorize -ftree-loop-vectorize'
		flags_o02='-march=native'

	flags=${flg_o01}' '${flg_o02}

# execution
	
	# remove modules and object codes & results.dat existing
	rm -f *.mod *.o
	rm -f result_y1.dat result_y2.dat result_energies.dat

	gfortran ${flags} -o ${ob_cod_name} ${pth_mod} ${f90_cod_name}
	
#	perf stat -e cycles,instructions,cache-references,cache-misses -r 5 ./${ob_cod_name}
	
	./${ob_cod_name}

	# remove modules and object codes
	rm -f *.mod *.o

# Notes:

# compilation:
# 	1) chmod +x script.sh
#	2) ./script.sh
