#!/bin/bash

# module paths
	pth_mod_01='../../../modules/module_presition.f90'
	pth_mod_02='../../../modules/module_functions_1D.f90'
	pth_mod_03='../../../modules/module_functions_2D.f90'
	pth_mod_04='../../../modules/module_EDO_segundo_orden.f90'
	pth_mod_05='../../../modules/module_numerical_error.f90'
	pth_mod_06='../../../modules/module_fourier_transform.f90'
	pth_mod_07='../../../modules/module_double_pendulum.f90'
	pth_mod_08='../../../modules/module_EDO_segundo_orden_flip.f90'

#	pth_mod=''${pth_mod_01}' '${pth_mod_03}' '${pth_mod_07}' '${pth_mod_04}' '${pth_mod_06}'' 	# RUN 2
	pth_mod=''${pth_mod_01}' '${pth_mod_07}' '${pth_mod_08}'' 					# RUN 3

# object code name
#	ob_cod_name='double_pendulum_spectrum.o'	# RUN 2
	ob_cod_name='heyl_double_pendulum.o'		# RUN 3

# fortran code name
#	f90_cod_name='double_pendulum_spectrum.f90' # RUN 2
	f90_cod_name='heyl_double_pendulum.f90' 	# RUN 3

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
#	lib01='-lfftw3' # RUN 2

# execution
	
	# remove modules and object codes & results.dat existing
	rm -f *.mod *.o
	
#	ld_path='-L/usr/local/lib/' # RUN 2
	
#	gfortran ${flags} -o ${ob_cod_name} ${pth_mod} ${f90_cod_name} ${ld_path} ${lib01} 	# RUN 2
	gfortran ${flags} -o ${ob_cod_name} ${pth_mod} ${f90_cod_name} 						# RUN 3
	
	./${ob_cod_name}

	# remove modules and object codes
	rm -f *.mod *.o

	./script_graph_gnuplot.sh
	
	rm -f result_flips.dat

# Notes:

# compilation:
# 	1) chmod +x script.sh
#	2) ./script.sh
