# ============================================================
# Especifico la regla que se ejuta por defecto
# ============================================================
#all: $(bin_name_01).o

# ============================================================
# COMPILERS
# ============================================================
# Fortran compiler
fcomp_01 = gfortran
# used fortran compiler
FC = $(fcomp_01)

# ============================================================
# COMPILER FLAGS
# ============================================================
# warning flags
flg_war01 = -Wall -Wextra -Wconversion -pedantic -ffpe-trap=overflow
# debugging flags
flg_deb01 = -ffpe-trap=zero -ffpe-trap=overflow -ffpe-trap=underflow 
flg_deb02 = -g -fcheck=all -fbacktrace -fbounds-check
# optimization flags
flg_opt01 = -o3 -ftree-vectorize -ftree-loop-vectorize
flg_opt02 = -march=native
# used compiler flags
FFLAGS = $(flg_opt01) $(flg_opt02) $(flg_war01) $(flg_deb01) $(flg_deb02)

# library flags
flg_lib01 = -lfftw3
# used library flags
LDFLAGS = $(flg_lib01)

# ============================================================
# BINARIES COMPILATION
# ============================================================
# binaries names
bin_name_01 = molecular_dynamic_lennard_jones_01
bin_name_02 = molecular_dynamic_lennard_jones_02
bin_name_03 = molecular_dynamic_lennard_jones_03
bin_name_04 = molecular_dynamic_lennard_jones_04
bin_name_05 = molecular_dynamic_lennard_jones_05
bin_name_06 = molecular_dynamic_lennard_jones_06
bin_name_07 = molecular_dynamic_lennard_jones_01_linked_list

# object dependencies
dependencies_01= module_precision.o module_mt19937.o module_md_lennard_jones.o

# binaries rules
# molecular_dynamic_lennard_jones_01 (inciso a)
$(bin_name_01).o: $(dependencies_01)
	$(FC) $(FFLAGS) -o $(bin_name_01).o $(dependencies_01) $(bin_name_01).f90

# molecular_dynamic_lennard_jones_02 (iniciso b)
$(bin_name_02).o: $(dependencies_01)
	$(FC) $(FFLAGS) -o $(bin_name_02).o $(dependencies_01) $(bin_name_02).f90

# molecular_dynamic_lennard_jones_03 (iniciso c)
$(bin_name_03).o: $(dependencies_01)
	$(FC) $(FFLAGS) -o $(bin_name_03).o $(dependencies_01) $(bin_name_03).f90

# molecular_dynamic_lennard_jones_04 (iniciso e)
$(bin_name_04).o: $(dependencies_01)
	$(FC) $(FFLAGS) -o $(bin_name_04).o $(dependencies_01) $(bin_name_04).f90

# molecular_dynamic_lennard_jones_05 (iniciso f)
$(bin_name_05).o: $(dependencies_01)
	$(FC) $(FFLAGS) -o $(bin_name_05).o $(dependencies_01) $(bin_name_05).f90

# molecular_dynamic_lennard_jones_06 (iniciso g)
$(bin_name_06).o: $(dependencies_01)
	$(FC) $(FFLAGS) -o $(bin_name_06).o $(dependencies_01) $(bin_name_06).f90

# molecular_dynamic_lennard_jones_01_linked_list
$(bin_name_07).o: $(dependencies_01)
	$(FC) $(FFLAGS) -o $(bin_name_07).o $(dependencies_01) $(bin_name_07).f90

# ============================================================
# MODULES COMPILATION
# ============================================================
# modules paths
mod_path = ../../../modules/
# modules names
mod_name_01 = module_precision
mod_name_02 = module_functions_1D
mod_name_03 = module_functions_2D
mod_name_04 = module_EDO_segundo_orden
mod_name_05 = module_numerical_error
mod_name_06 = module_fourier_transform
mod_name_07 = module_double_pendulum
mod_name_08 = module_EDO_segundo_orden_flip
mod_name_19 = module_pullen_edmonds
mod_name_10 = module_EDO_segundo_orden_poincare
mod_name_11 = module_tridiag_matrix
mod_name_12 = module_random_generator
mod_name_13 = module_mzran
mod_name_14 = module_mt19937
mod_name_15 = module_2D_ferromagnetic_ising
mod_name_16 = module_md_lennard_jones

# modules recipes
$(mod_name_01).o:
	$(FC) -c $(mod_path)$(mod_name_01).f90
$(mod_name_12).o: $(mod_name_01).o
	$(FC) -c $(mod_path)$(mod_name_12).f90
$(mod_name_13).o: $(mod_name_01).o
	$(FC) -c $(mod_path)$(mod_name_13).f90
$(mod_name_14).o:
	$(FC) -c $(mod_path)$(mod_name_14).f90
$(mod_name_15).o: $(mod_name_01).o $(mod_name_14).o
	$(FC) -c $(mod_path)$(mod_name_15).f90
$(mod_name_16).o: $(mod_name_01).o $(mod_name_14).o
	$(FC) -c $(mod_path)$(mod_name_16).f90	
	
# ============================================================
# PONY rules
# ============================================================
clean:
	rm -f *.o *.mod

# References
# https://www.gnu.org/software/make/manual/make.html
# https://faculty.washington.edu/rjl/classes/am583s2013/notes/gfortran_flags.html