program estimated_size
    use module_precision
    implicit none
    real(dp) :: valor1=4._dp*atan(1.0_dp),valor2=4._dp*atan(1.0_dp)
    integer(sp) :: i,istat
    open(10,file='estimated_size.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'Input/Output file. istat10 = ',istat
    20 format(E9.2,x,E9.2)
    do i=1,200; write(10,20) valor1,valor2; end do
end program estimated_size

! gfortran -o estimated_size.o module_precision.f90 estimated_size.f90 && ./estimated_size.o && rm *.mod *.o