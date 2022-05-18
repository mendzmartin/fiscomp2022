! make clean && make test_integer_rnd.o && ./test_integer_rnd.o
program test_integer_rnd
    use module_precision;use module_random_generator
    implicit none
    integer(sp), parameter :: m=1,n=10
    integer(sp)            :: seed,seed_val(8)
    integer(sp)            :: i,j,i_rnd,j_rnd
    real(dp)               :: nrand
    ! generate seed
    call date_and_time(values=seed_val)
    seed=seed_val(8)*seed_val(7)*seed_val(6)+seed_val(5) 
    ! floor(x,type) returns the greatest integer less than or equal to X
    do i=1,n;do j=1,n
        ! We want to choose one from m to n integers
        nrand=ran2(seed);i_rnd=m+floor((n+1-m)*nrand,sp) ! rows
        nrand=ran2(seed);j_rnd=m+floor((n+1-m)*nrand,sp) ! columns
        write(*,*) i_rnd,j_rnd
    end do;end do
end program test_integer_rnd