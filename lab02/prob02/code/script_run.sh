#!/bin/bash

# module paths
	pth_mod_01='../../../modules/module_precision.f90'
	pth_mod_02='../../../modules/module_functions_1D.f90'
	pth_mod_03='../../../modules/module_functions_2D.f90'
	pth_mod_04='../../../modules/module_EDO_segundo_orden.f90'
	pth_mod_05='../../../modules/module_numerical_error.f90'
	pth_mod_06='../../../modules/module_fourier_transform.f90'
	pth_mod_07='../../../modules/module_double_pendulum.f90'
	pth_mod_08='../../../modules/module_EDO_segundo_orden_flip.f90'
	pth_mod_09='../../../modules/module_pullen_edmonds.f90'
	pth_mod_10='../../../modules/module_EDO_segundo_orden_poincare.f90'
	pth_mod_11='../../../modules/module_tridiag_matrix.f90'

	pth_mod=''${pth_mod_01}'' 				  # RUN 1,4
#	pth_mod=''${pth_mod_01}' '${pth_mod_11}'' # RUN 2,3,5,6

# object code name
#	ob_cod_name='heateq_explicit_von_neumann.o' 	  #RUN 1
#	ob_cod_name='heateq_implicit_von_neumann.o'       #RUN 2
#	ob_cod_name='heateq_crank_nicolson_von_neumann.o' #RUN 3
	ob_cod_name='heateq_comparison_01.o' 		      #RUN 4
#	ob_cod_name='heateq_comparison_02.o'              #RUN 5
#	ob_cod_name='heateq_comparison_03.o' 			  #RUN 6


# fortran code name
#	f90_cod_name='heateq_explicit_von_neumann.f90' 		 #RUN 1
#	f90_cod_name='heateq_implicit_von_neumann.f90' 		 #RUN 2
#	f90_cod_name='heateq_crank_nicolson_von_neumann.f90' #RUN 3
	f90_cod_name='heateq_comparison_01.f90' 		     #RUN 4
#	f90_cod_name='heateq_comparison_02.f90' 			 #RUN 5
#	f90_cod_name='heateq_comparison_03.f90'				 #RUN 6

#flags compiles
#{see: https://faculty.washington.edu/rjl/classes/am583s2013/notes/gfortran_flags.html}
	## warning flags
		flg_w01='-Wall -Wextra -Wconversion -pedantic'
	## debugging flags
		flg_d01='-ffpe-trap=zero -ffpe-trap=overflow -ffpe-trap=underflow' 
		flg_d02='-g -fbacktrace -fbounds-check'
	## optimization flags
		flg_o01='-o3 -ftree-vectorize -ftree-loop-vectorize'
		flg_o02='-march=native'

	flags=${flg_o01}' '${flg_o02}' '${flg_w01}' '${flg_d02}
	
# libraries
#	lib01='-lfftw3'
# execution
	
	# remove modules and object codes & results.dat existing
	rm -f *.mod *.o
	
#	ld_path='-L/usr/local/lib/'
	
#	gfortran ${flags} -o ${ob_cod_name} ${pth_mod} ${f90_cod_name} ${ld_path} ${lib01}
	gfortran ${flags} -o ${ob_cod_name} ${pth_mod} ${f90_cod_name} 						# RUN 1,2,3,4
	
	./${ob_cod_name}

	# remove modules and object codes
	rm -f *.mod *.o

# Notes:
# perf stat -e cycles,instructions,cache-references,cache-misses -r 20 ./script_run.sh

# compilation:
# 	1) chmod +x script.sh
#	2) ./script.sh
