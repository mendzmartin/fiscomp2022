! module of monte carlo to lennard jones potential
! gfortran -c module_precision.f90 module_mt19937.f90 module_md_lennard_jones.f90
module module_md_lennard_jones
    use module_precision;use module_mt19937, only: sgrnd,grnd
    implicit none
    contains

    subroutine init_molecular_dynamic(n_p,x_vector,y_vector,z_vector,U_adim,T_adim,&
                                      r_cutoff,density)
        implicit none
        integer(sp), intent(in)    :: n_p                                       ! cantidad de partículas
        real(dp),    intent(inout) :: x_vector(n_p),y_vector(n_p),z_vector(n_p) ! vectores posicion
        real(dp),    intent(inout) :: U_adim                                    ! Energía interna
        real(dp),    intent(in)    :: T_adim                                    ! Temperatura de equilibrio
        real(dp),    intent(in)    :: r_cutoff ! distancia de truncado
        real(dp),    intent(in)    :: density
        
        real(dp)    :: U_adim_new,U_adim_old,deltaU_adim ! energías con/sin despl. y variación de energía
        real(dp)    :: x_old,y_old,z_old                 ! posiciones sin deplazar
        real(dp)    :: delta_x,delta_y,delta_z           ! desplazamientos
        real(dp)    :: nrand                             ! numero random
        integer(sp) :: seed,seed_val(8),MC_index,index

        call date_and_time(values=seed_val)
        seed=seed_val(8)*seed_val(7)*seed_val(6)+seed_val(5);call sgrnd(seed)
        delta_x=0.05_dp;delta_y=0.05_dp;delta_z=0.05_dp ! ¿CÓMO DEFINO ESTOS VALORES?

        ! hago un paso de Monte Carlo (MC_step)
        do MC_index=1,n_p
            U_adim_old=U_adim ! energía sin desplazar
            ! random integer from 1 to n_p (elijo partic al azar y desplazo)
            nrand=real(grnd(),dp);index=1_sp+floor(n_p*nrand,sp)
            x_old=x_vector(index);y_old=y_vector(index);z_old=z_vector(index)
            nrand=real(grnd(),dp);x_vector(index)=x_old+(nrand-0.5_dp)*delta_x
            nrand=real(grnd(),dp);y_vector(index)=y_old+(nrand-0.5_dp)*delta_y
            nrand=real(grnd(),dp);z_vector(index)=z_old+(nrand-0.5_dp)*delta_z
            ! corrección de posiciones
            call position_correction(n_p,density,x_vector(index),y_vector(index),z_vector(index))
            ! calculo energía desplazada
            U_adim_new=u_lj_total(n_p,x_vector,y_vector,z_vector,r_cutoff)
            deltaU_adim=(U_adim_new-U_adim_old) ! variación de energía interna
            cond1:  if (deltaU_adim<=0._dp) then;U_adim=U_adim_new;exit cond1
                    else
                        if (T_adim==0._dp) exit cond1
                        nrand=real(grnd(),dp)
                        if ((deltaU_adim*(1._dp/T_adim))<abs(log(tiny(1._dp)))) then
                            if (exp(-deltaU_adim*(1._dp/T_adim))>=nrand) U_adim=U_adim_new
                            exit cond1
                        else
                            if (nrand==0._dp) then;U_adim=U_adim_new;exit cond1;end if
                            U_adim=U_adim_old
                            x_vector(index)=x_old;y_vector(index)=y_old;z_vector(index)=z_old
                        end if
                    end if cond1
        end do
    end subroutine

    ! corrección de las posiciones (PBC)
    subroutine position_correction(n_p,density,x,y,z)
        integer(sp), intent(in)    :: n_p       ! numero total de partículas
        real(dp),    intent(in)    :: density   ! densidad de partículas
        real(dp),    intent(inout) :: x,y,z     ! posiciones
        real(dp)                   :: L         ! logitud macroscópica por dimensión
        L=(real(n_p,dp)*(1._dp/density))**(1._dp/3._dp)
        if (x<0._dp) x=x+L;if (y<0._dp) y=y+L;if (z<0._dp) z=z+L
        if (x>L) x=x-L;if (y>L) y=y-L;if (z>L) z=z-L
    end subroutine position_correction

    ! Subroutine to set up fcc lattice
    subroutine initial_particle_configuration_fcc(n_p,density,x_vector,y_vector,z_vector)
        integer(sp), intent(in)    :: n_p                ! numero total de partículas
        real(dp),    intent(in)    :: density            ! densidad de partículas
        real(dp),    intent(inout) :: x_vector(n_p),&    ! coordenadas de los vectores posición
                                      y_vector(n_p),&
                                      z_vector(n_p)
        integer(sp)                :: i,j,k,index2,index ! loop index
        real(dp)                   :: L                  ! logitud macroscópica por dimensión
        real(dp)                   :: a                  ! parámetro de red
        real(dp)                   :: n_unitcells        ! numero de celdas unidad por dimension
        real(dp),   allocatable    :: aux_matrix(:,:)    ! matriz auxiliar de indices (FCC)
        L=(real(n_p,dp)*(1._dp/density))**(1._dp/3._dp)
        n_unitcells=anint((real(n_p,dp)*0.25_dp)**(1._dp/3._dp),dp)
        a=L*(1._dp/n_unitcells)!a=(4*(1._dp/density))**(1._dp/3._dp)
        ! cargamos datos en matrix auxiliar (specific for FCC structure)
        allocate(aux_matrix(4,3));aux_matrix(:,:)=a*0.5_dp
        aux_matrix(1,:)=0.0_dp;aux_matrix(2,3)=0.0_dp
        aux_matrix(3,2)=0.0_dp;aux_matrix(4,1)=0.0_dp
        ! cargamos vectores de coordenadas
        index=0
        do i=1,int(n_unitcells,sp);do j=1,int(n_unitcells,sp);do k=1,int(n_unitcells,sp);do index2=1,4
            index=index+1
            x_vector(index)=aux_matrix(index2,1)+real(i-1,dp)*a
            y_vector(index)=aux_matrix(index2,2)+real(j-1,dp)*a
            z_vector(index)=aux_matrix(index2,3)+real(k-1,dp)*a
        end do;end do;end do;end do
        deallocate(aux_matrix)
    end subroutine initial_particle_configuration_fcc

    ! compute individual lennard jones potential (simple truncation)
    function u_lj_individual(x1,y1,z1,x2,y2,z2,r_cutoff)
        real(dp), intent(in) :: x1,y1,z1,x2,y2,z2 ! coordenadas del par de partículas
        real(dp), intent(in) :: r_cutoff
        real(dp)             :: factor1,factor2   ! factores posiciones
        real(dp)             :: r12               ! distancia adimensional entre pares de particulas
        real(dp)             :: u_lj_individual   ! adimensional lennard jones potential
        integer(sp)          :: i
        r12=abs(sqrt(x1*x1+y1*y1+z1*z1)-sqrt(x2*x2+y2*y2+z2*z2))
        if (r12/=0._dp) then
            if (r12<=r_cutoff) then
                factor1=1._dp;factor2=1._dp
                do i=1,12;factor1=factor1*(1._dp/r12);end do ! (r12)^12
                do i=1,6;factor2=factor2*(1._dp/r12);end do  ! (r12)^6
                u_lj_individual=4*(factor1-factor2)
            else;u_lj_individual=0._dp
            end if
        else
            u_lj_individual=0._dp
        end if
    end function u_lj_individual
    ! compute total lennard jones potential
    function u_lj_total(n_p,x_vector,y_vector,z_vector,r_cutoff)
        integer(sp), intent(in) :: n_p
        real(dp),    intent(in) :: x_vector(n_p),y_vector(n_p),z_vector(n_p)
        real(dp),    intent(in) :: r_cutoff
        real(dp)                :: u_lj_total
        integer(sp)             :: i,j
        u_lj_total=0._dp
        do j=1,n_p
            do i=1,j-1
                u_lj_total=u_lj_total+u_lj_individual(x_vector(i),y_vector(i),z_vector(i),&
                                                  x_vector(j),y_vector(j),z_vector(j),r_cutoff)
            end do
        end do
    end function u_lj_total

end module module_md_lennard_jones