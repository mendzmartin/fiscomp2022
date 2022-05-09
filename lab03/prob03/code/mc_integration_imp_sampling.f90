! problema 03.b
! ojo la distribución debe normalizarse segun estos valores
! este programa sólo vale para cuando x_start=0; x_end=1
program mc_integration_imp_sampling
    use module_precision;use module_random_generator
    implicit none
    
    integer(sp), parameter :: pot=3_sp   ! potencia (debe ser mayor a -2)
    integer(sp), parameter :: potk=3_sp  ! usar valores 2 y 3
    integer(sp)            :: n ! cantidad de evaluaciones de la función
    real(dp),    parameter :: x_start=0._dp,x_end=1._dp
    real(dp),    parameter :: exact_integ=1._dp/(real(pot,dp)+1._dp)
    integer(sp)            :: seed,seed_val(8),i,j,k,istat
    real(dp)               :: nrand,x_rand,f_rand,g_rand,integ
    real(dp)               :: factor

    call date_and_time(values=seed_val)
    seed=seed_val(8)*seed_val(7)*seed_val(6)+seed_val(5)

    integ=0._dp
    !open(10,file='../results/result_P03b_01.dat',status='replace',action='write',iostat=istat) ! p/ potk=2
    open(10,file='../results/result_P03b_02.dat',status='replace',action='write',iostat=istat)  ! p/ potk=3
    if (istat /= 0_sp) write(*,*) 'istat_error=', istat
    20 format(2(A12,x),A12);21 format(I12,x,E12.4,x,E12.4)
    write(10,20) 'n','Iaprox','E_rel'
    do i=1,1E+04
            n=10_sp*i
        do j=1,n
            nrand=ran0(seed)
            ! x^k distribution (x \in {1,Infinity})
            x_rand=nrand**(1_dp/real(potk+1_dp,dp))
            f_rand=1._dp;do k=1,pot;f_rand=f_rand*x_rand;end do  ! x_rand**pot
            factor=1._dp;do k=1,potk;factor=factor*x_rand;end do ! x_rand**potk
            g_rand=(real(potk,dp)+1_dp)*factor
            integ=integ+f_rand*(1._dp/g_rand)
        end do
        integ=(x_end-x_start)*(1._dp/real(n))*integ
        write(10,21) n,integ,abs((exact_integ-integ)*(1._dp/exact_integ))
    end do
    close(10)

    ! controlamos valores de la integral en el último paso
    write(*,'(A10,E12.4)') 'Iaprox=',integ
    write(*,'(A10,E12.4)') 'Iexact=',1._dp/(real(pot,dp)+1._dp)

end program mc_integration_imp_sampling