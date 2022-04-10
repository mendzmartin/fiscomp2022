#!/bin/bash

# module paths
	pth_mod_01='../../../modules/module_presition.f90'
	pth_mod_02='../../../modules/module_functions_1D.f90'
	pth_mod_03='../../../modules/module_functions_2D.f90'
	pth_mod_04='../../../modules/module_EDO_segundo_orden.f90'
	pth_mod_05='../../../modules/module_numerical_error.f90'
	pth_mod_06='../../../modules/module_fourier_transform.f90'
	pth_mod_07='../../../modules/module_double_pendulum.f90'

	pth_mod=''${pth_mod_01}' '${pth_mod_03}' '${pth_mod_04}''

# object code name
	#ob_cod_name='logistic_map.o'
	ob_cod_name='hola.o'

# fortran code name
	#f90_cod_name='logistic_map.f90'
	f90_cod_name='hola.f90'

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

	flags=${flg_o01}' '${flg_o02}' '${flg_w01}' '${flg_d02} 
	
# libraries
	#lib01='-lfftw3'

# execution
	
	# remove modules and object codes & results.dat existing
	rm -f *.mod *.o
	
	#export LD_LIBRARY_PATH= '-L../../../libraries/'
	#ld_path='-L/usr/local/lib/'
	
	#gfortran ${flags} -o ${ob_cod_name} ${pth_mod} ${f90_cod_name} ${ld_path} ${lib01}
	gfortran ${flags} -o ${ob_cod_name} ${pth_mod} ${f90_cod_name}
	
	./${ob_cod_name}

	# remove modules and object codes
	rm -f *.mod *.o

# Notes:

# compilation:
# 	1) chmod +x script.sh
#	2) ./script.sh
