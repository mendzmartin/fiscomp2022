program test_module_mt19937
    ! this main() outputs first 1000 generated numbers
    use module_mt19937
    implicit none
    
    integer, parameter      :: no=1000
    real(8), dimension(0:7) :: r
    integer j,k
    ! real(8) grnd
    ! call sgrnd(435754)
    ! write(*,*) 'mt=',grnd()
    any nonzero integer can be used as a seed
    do j=0,no-1
        r(mod(j,8))=grnd()
        if(mod(j,8).eq.7) then
        write(*,'(8(f9.6,'' ''))') (r(k),k=0,7)
        else if(j.eq.no-1) then
        write(*,'(8(f9.6,'' ''))') (r(k),k=0,mod(no-1,8))
        endif
    enddo
    
end program test_module_mt19937
! run: gfortran -Wall -o test_module_mt19937.o ../module_mt19937.f90 test_module_mt19937.f90
!      ./test_module_mt19937.o && rm *.mod *.o