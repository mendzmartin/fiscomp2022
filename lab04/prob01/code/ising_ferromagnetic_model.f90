program ising_ferromagnetic_model
    use module_precision;use module_random_generator
    implicit none

    integer(sp), parameter   :: n=1000_sp ! sitios de red por dimension
    real(dp), parameter      :: J_int=1._sp ! J_int>0 => ferromagnetic
    integer(sp), allocatable :: spin_matrix(:,:)
    real(dp) :: T_adim ! temperatura adimensional
    integer(sp) :: seed,seed_val(8),i,j
    real(dp)    :: nrand


    ! genero configuracion inicial
    call date_and_time(values=seed_val)
    seed=seed_val(8)*seed_val(7)*seed_val(6)+seed_val(5)

    allocate(spin_matrix(n,n))
    do j=1,n ! recorro columnas
        do i=1,n ! recorro filas
            nrand=ran2(seed)
            if (nrand>=0._dp.and.nrand<0.5_dp) then;spin_matrix(i,j)=1_sp
            else; spin_matrix(i,j)=0_sp;end if
        end do
    end do

end program ising_ferromagnetic_model