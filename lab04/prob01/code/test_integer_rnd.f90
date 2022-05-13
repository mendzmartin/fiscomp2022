program test_integer_rnd
    use module_precision
    use module_random_generator
    implicit none

    integer(sp), parameter :: m=1,n=10
    integer(sp)            :: seed,seed_val(8),i,j,i_rnd,j_rnd
    real(dp)               :: nrand

    call date_and_time(values=seed_val)
    seed=seed_val(8)*seed_val(7)*seed_val(6)+seed_val(5)
    
    ! floor(x,type) returns the greatest integer less than or equal to X

    do i=1,n
        do j=1,n
            nrand=ran2(seed)
            i_rnd=m+floor((n+1-m)*nrand,sp) ! We want to choose one from m to n integers
            nrand=ran2(seed)
            j_rnd=m+floor((n+1-m)*nrand,sp)
            write(*,*) i_rnd,j_rnd
        end do
    end do

end program test_integer_rnd