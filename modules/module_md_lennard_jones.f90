! module of molecular dynamic to lennard jones potential
module module_md_lennard_jones
    use module_precision;use module_mt19937, only: sgrnd,grnd
    implicit none
    contains

    subroutine init_molecular_dynamic(n_p,x_vector,y_vector,z_vector,U_adim)
        implicit none
        integer(sp), intent(in)    :: n_p            ! Tipo de MC relaxation y dimension
        integer(sp), intent(inout) :: x_vector(n_p),y_vector(n_p),z_vector(n_p)   ! Matriz de spins
        real(dp),    intent(inout) :: U_adim                    ! Energía interna
        
        real(dp)    :: nrand
        real(dp)    :: U_adim_new,U_adim_old
        integer(sp) :: seed,seed_val(8),i,j,k,index

        call date_and_time(values=seed_val)
        seed=seed_val(8)*seed_val(7)*seed_val(6)+seed_val(5);call sgrnd(seed)

        ! hago un paso de Monte Carlo (MC_step)
        do k=1,n_p
            U_adim_old=u_lj_total(x_vector,y_vector,z_vector)
            nrand=ran2(seed);index=1_sp+floor(n_p*nrand,sp) ! random integer from 1 to n_p (elijo partic al azar)

            ! cond1:if (deltaU_adim<=0._dp) then
            !     U_adim=U_adim+deltaU_adim
            !     aux_matrix_pbc(i,j)=-aux_matrix_pbc(i,j)
            !     call pbc_ij(aux_matrix_pbc,n,i,j) ! actualizamos condiciones de borde
            ! else
            !     if (T_adim==0._dp) exit cond1
            !     nrand=real(grnd(),dp)
            !     if (deltaU_adim*(1._dp/T_adim)<abs(log(tiny(1._dp)))) then
            !         if (exp(-deltaU_adim*(1._dp/T_adim))>=nrand) then
            !             U_adim=U_adim+deltaU_adim
            !             aux_matrix_pbc(i,j)=-aux_matrix_pbc(i,j)
            !             call pbc_ij(aux_matrix_pbc,n,i,j) ! actualizamos condiciones de borde
            !         end if
            !     else
            !         if (nrand==0._dp) then
            !             U_adim=U_adim+deltaU_adim
            !             aux_matrix_pbc(i,j)=-aux_matrix_pbc(i,j)
            !             call pbc_ij(aux_matrix_pbc,n,i,j) ! actualizamos condiciones de borde
            !         end if
            !     end if
            ! end if cond1
        end do
    end subroutine

    ! Subroutine to set up fcc lattice
    subroutine initial_particle_configuration_fcc(n_p,density,x_coord_vector,y_coord_vector,z_coord_vector)
        integer(sp), intent(in)    :: n_p                   ! numero total de partículas
        real(dp),    intent(in)    :: density               ! densidad de partículas
        real(dp),    intent(inout) :: x_coord_vector(n_p),& ! coordenadas de los vectores posición
                                      y_coord_vector(n_p),&
                                      z_coord_vector(n_p)
        integer(sp)           :: i,j,k,index2,index ! loop index
        real(dp)              :: L                  ! logitud macroscópica por dimensión
        real(dp)              :: a                  ! parámetro de red
        real(dp)              :: n_unitcells        ! numero de celdas unidad por dimension
        real(dp), allocatable :: aux_matrix(:,:)    ! matriz auxiliar de indices (FCC)
        L=(n_p*(1._dp/density))**(1._dp/3._dp)   
        n_unitcells=anint((n_p*0.25_dp)**(1._dp/3._dp),dp) 
        a=L*(1._dp/n_unitcells)!a=(4*(1._dp/density))**(1._dp/3._dp)
        ! cargamos datos en matrix auxiliar (FCC)
        allocate(aux_matrix(4,3));aux_matrix(:,:)=a*0.5_dp
        aux_matrix(1,:)=0.0_dp;aux_matrix(2,3)=0.0_dp
        aux_matrix(3,2)=0.0_dp;aux_matrix(4,1)=0.0_dp
        ! cargamos vectores de coordenadas
        index=0
        do i=1,int(n_unitcells,sp);do j=1,int(n_unitcells,sp);do k=1,int(n_unitcells,sp);do index2=1,4
            index=index+1
            x_coord_vector(index)=aux_matrix(index2,1)+real(i-1,dp)*a
            y_coord_vector(index)=aux_matrix(index2,2)+real(j-1,dp)*a
            z_coord_vector(index)=aux_matrix(index2,3)+real(k-1,dp)*a
        end do;end do;end do;end do
        deallocate(aux_matrix)
    end subroutine initial_particle_configuration_fcc

    ! compute individual lennard jones potential
    function u_lj_individual(x1,y1,z1,x2,y2,z2)
        real(dp), intent(in) :: x1,y1,z1,x2,y2,z2 ! coordenadas del par de partículas
        real(dp)             :: factor1,factor2   ! factores posiciones
        real(dp)             :: r12               ! distancia adimensional entre pares de particulas
        real(dp)             :: u_lj              ! adimensional lennard jones potential
        integer(sp)          :: i
        r12=abs(sqrt(x1*x1+y1*y1+z1*z1)-sqrt(x2*x2+y2*y2+z2*z2))
        factor1=1._dp;factor2=1._dp
        do i=1,12;factor1=factor1*(1._dp/r12);end do ! (r12)^12
        do i=1,6;factor2=factor2*(1._dp/r12);end do  ! (r12)^6
        u_lj=4*(factor1-factor2)
    end function u_lj_individual
    ! compute total lennard jones potential
    function u_lj_total(n_p,x_vector,y_vector,z_vector)
        integer(sp), intent(in) :: n_p
        real(dp),    intent(in) :: x_vector(np),y_vector(np),z_vector(np)
        real(dp)                :: u_lj_total
        integer(sp)             :: i
        u_lj_total=0._dp
        do i=1,n_p-1
            u_lj_total=u_lj_total+u_lj_individual(x_vector(i),y_vector(i),z_vector(i),&
                                                  x_vector(i+1),y_vector(i+1),z_vector(i+1))
        end do
    end function u_lj_total



end module module_md_lennard_jones