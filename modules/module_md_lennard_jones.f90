! module of molecular dynamic to lennard jones potential
! gfortran -c module_precision.f90 module_mt19937.f90 module_md_lennard_jones.f90
module module_md_lennard_jones
    use module_precision;use module_mt19937, only: sgrnd,grnd
    implicit none
    contains

    ! FUNCIONES
    ! compute individual lennard jones potential (simple truncation)
    function u_lj_individual(n_p,x1,y1,z1,x2,y2,z2,r_cutoff,density)
        integer(sp), intent(in) :: n_p
        real(dp), intent(in) :: x1,y1,z1,x2,y2,z2   ! coordenadas del par de partículas
        real(dp), intent(in) :: r_cutoff,density
        real(dp)             :: r12_pow06,r12_pow12 ! potencias de la distancia relativa
        real(dp)             :: r12                 ! distancia adimensional entre pares de particulas
        real(dp)             :: u_lj_individual     ! adimensional lennard jones potential
        integer(sp)          :: i
        ! calculamos distancia relativa corregida según PBC
        r12=rel_pos_correction(x1,y1,z1,x2,y2,z2,n_p,density)
        if (r12/=0._dp) then
            if (r12<=r_cutoff) then
                r12_pow06=1._dp
                do i=1,6;r12_pow06=r12_pow06*(1._dp/r12);end do  ! (r12)^6
                r12_pow12=r12_pow06*r12_pow06 ! (r12)^12
                u_lj_individual=4._dp*(r12_pow12-r12_pow06)
            else;u_lj_individual=0._dp
            end if
        else
            write(*,*) 'r12=0'
            u_lj_individual=0._dp
        end if
    end function u_lj_individual

    function pressure(n_p,density,mass,r_cutoff,x_vector,y_vector,z_vector,&
        vx_vector,vy_vector,vz_vector)
        integer(sp), intent(in) :: n_p
        real(dp),    intent(in) :: r_cutoff,density,mass
        real(dp),    intent(in) :: x_vector(n_p),y_vector(n_p),z_vector(n_p)
        real(dp),    intent(in) :: vx_vector(n_p),vy_vector(n_p),vz_vector(n_p)
        real(dp)                :: pressure,rij,result,T_adim
        integer(sp) :: i,j
        result=0._dp
        do j=1,n_p;do i=1,j-1
            rij=rel_pos_correction(x_vector(i),y_vector(i),z_vector(i),&
            x_vector(j),y_vector(j),z_vector(j),n_p,density)
            result=result+f_lj_individual(x_vector(i),y_vector(i),z_vector(i),&
            x_vector(j),y_vector(j),z_vector(j),r_cutoff,n_p,density)*rij
        end do;end do
        T_adim=temperature(n_p,mass,vx_vector,vy_vector,vz_vector)
        pressure=density*(T_adim+(1._dp/(3._dp*real(n_p,dp)))*result)
    end function pressure

    ! compute total lennard jones potential
    function u_lj_total(n_p,x_vector,y_vector,z_vector,r_cutoff,density)
        integer(sp), intent(in) :: n_p
        real(dp),    intent(in) :: x_vector(n_p),y_vector(n_p),z_vector(n_p)
        real(dp),    intent(in) :: r_cutoff,density
        real(dp)                :: u_lj_total
        integer(sp)             :: i,j
        u_lj_total=0._dp
        do j=1,n_p
            do i=1,j-1
                u_lj_total=u_lj_total+u_lj_individual(n_p,x_vector(i),y_vector(i),z_vector(i),&
                x_vector(j),y_vector(j),z_vector(j),r_cutoff,density)
            end do
        end do
    end function u_lj_total

    function kinetic_ergy_total(n_p,vx_vector,vy_vector,vz_vector,mass)
        integer(sp), intent(in) :: n_p
        real(dp),    intent(in) :: vx_vector(n_p),vy_vector(n_p),vz_vector(n_p)
        real(dp),    intent(in) :: mass
        real(dp)                :: kinetic_ergy_total
        kinetic_ergy_total=0.5_dp*mass*(sum(vx_vector(:)*vx_vector(:))+&
        sum(vy_vector(:)*vy_vector(:))+sum(vz_vector(:)*vz_vector(:)))*(1._dp/real(n_p,dp))
    end function kinetic_ergy_total

    ! caclulo de la fuerza individual (par de partículas)
    function f_lj_individual(x1,y1,z1,x2,y2,z2,r_cutoff,n_p,density)
        real(dp),    intent(in) :: x1,y1,z1,x2,y2,z2 ! coordenadas del par de partículas
        integer(sp), intent(in) :: n_p ! cantidad total de partículas
        real(dp),    intent(in) :: r_cutoff,density
        real(dp)    :: r12_pow02,r12_pow06,r12_pow12   ! factores potencia
        real(dp)    :: r12               ! distancia adimensional entre pares de particulas
        real(dp)    :: f_lj_individual   ! adimensional individual lennard jones force
        integer(sp) :: i
        ! calculamos distancia relativa corregida según PBC
        r12=rel_pos_correction(x1,y1,z1,x2,y2,z2,n_p,density)
        r12_pow02=r12*r12
        if (r12/=0._dp) then
            if (r12<=r_cutoff) then
                r12_pow06=1._dp
                do i=1,3;r12_pow06=r12_pow06*(1._dp/r12_pow02);end do ! (r12)^6
                r12_pow12=r12_pow06*r12_pow06 ! (r12)^12
                f_lj_individual=48._dp*(1._dp/r12_pow02)*(r12_pow12-0.5_dp*r12_pow06)
            else;f_lj_individual=0._dp
            end if
        else
            write(*,*) 'r12=0'
            f_lj_individual=0._dp
        end if
    end function f_lj_individual

    ! calculo de la componente xi de la fuerza total
    function f_lj_total(x_vector,y_vector,z_vector,r_cutoff,n_p,density)
        integer(sp), intent(in)    :: n_p
        real(dp),    intent(in)    :: x_vector(n_p),y_vector(n_p),z_vector(n_p)
        real(dp),    intent(in)    :: r_cutoff,density
        real(dp)                   :: f_lj_total
        integer(sp)                :: i,j
        f_lj_total=0._dp
        do j=1,n_p
            do i=1,j-1
                f_lj_total=f_lj_total+f_lj_individual(x_vector(i),y_vector(i),z_vector(i),&
                x_vector(j),y_vector(j),z_vector(j),r_cutoff,n_p,density)
            end do
        end do
    end function f_lj_total

    function temperature(n_p,mass,vx_vector,vy_vector,vz_vector)
        integer(sp), intent(in) :: n_p
        real(dp),    intent(in) :: vx_vector(n_p),vy_vector(n_p),vz_vector(n_p)
        real(dp),    intent(in) :: mass
        real(dp)                :: temperature
        temperature=mass*(sum(vx_vector(:)*vx_vector(:))+&
        sum(vy_vector(:)*vy_vector(:))+sum(vz_vector(:)*vz_vector(:)))*&
        (1._dp/(3._dp*real(n_p,dp)))
    end function 
    
    ! corrección de las posiciones relativas (PBC)
    function rel_pos_correction(x1,y1,z1,x2,y2,z2,n_p,density)
        integer(sp), intent(in)    :: n_p               ! numero total de partículas
        real(dp),    intent(in)    :: density           ! densidad de partículas
        real(dp),    intent(in)    :: x1,y1,z1,x2,y2,z2 ! coordenadas del par de partículas
        real(dp)                   :: rel_pos_correction ! posición relativa corregida (PBC)
        real(dp)                   :: L                 ! logitud macroscópica por dimensión
        real(dp)                   :: dx,dy,dz          ! diferencia de posiciones segun coordenadas
        L=(real(n_p,dp)*(1._dp/density))**(1._dp/3._dp)
        dx=abs(x2-x1);dy=abs(y2-y1);dz=abs(z2-z1)
        dx=dx+0.5_dp*L;dx=dx-L*anint(dx*(1._dp/L),dp);dx=abs(dx-0.5_dp*L)
        dy=dy+0.5_dp*L;dy=dy-L*anint(dy*(1._dp/L),dp);dy=abs(dy-0.5_dp*L)
        dz=dz+0.5_dp*L;dz=dz-L*anint(dz*(1._dp/L),dp);dz=abs(dz-0.5_dp*L)
        rel_pos_correction=sqrt(dx*dx+dy*dy+dz*dz)
    end function rel_pos_correction

    ! SUBRUTINAS
    subroutine rescaling_velocities(n_p,vx_vector,vy_vector,vz_vector,T_adim_ref,mass)
        integer(sp), intent(in)    :: n_p
        real(dp),    intent(inout) :: vx_vector(n_p),vy_vector(n_p),vz_vector(n_p)
        real(dp),    intent(in)    :: T_adim_ref,mass
        real(dp)                   :: T_adim
        T_adim=temperature(n_p,mass,vx_vector,vy_vector,vz_vector)
        vx_vector(:)=sqrt(T_adim_ref*(1._dp/T_adim))*vx_vector(:)
        vy_vector(:)=sqrt(T_adim_ref*(1._dp/T_adim))*vy_vector(:)
        vz_vector(:)=sqrt(T_adim_ref*(1._dp/T_adim))*vz_vector(:)
    end subroutine rescaling_velocities

    ! corrección de las posiciones (PBC)
    subroutine position_correction(n_p,density,x,y,z)
        integer(sp), intent(in)    :: n_p       ! numero total de partículas
        real(dp),    intent(in)    :: density   ! densidad de partículas
        real(dp),    intent(inout) :: x,y,z     ! posiciones
        real(dp)                   :: L         ! logitud macroscópica por dimensión
        L=(real(n_p,dp)*(1._dp/density))**(1._dp/3._dp)
        x=x+0.5_dp*L;x=x-L*anint(x*(1._dp/L),dp);x=x-0.5_dp*L
        y=y+0.5_dp*L;y=y-L*anint(y*(1._dp/L),dp);y=y-0.5_dp*L
        z=z+0.5_dp*L;z=z-L*anint(z*(1._dp/L),dp);z=z-0.5_dp*L
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
            ! CENTRAMOS LA CELDA EN EL RANGO [-L/2:L/2]
            x_vector(index)=(aux_matrix(index2,1)+real(i-1,dp)*a)-0.5_dp*L
            y_vector(index)=(aux_matrix(index2,2)+real(j-1,dp)*a)-0.5_dp*L
            z_vector(index)=(aux_matrix(index2,3)+real(k-1,dp)*a)-0.5_dp*L
        end do;end do;end do;end do
        deallocate(aux_matrix)
    end subroutine initial_particle_configuration_fcc

    ! SUBRUTINA DE INTEGRACIÓN DE ECUACIONES DE MOVIMIENTO
    subroutine velocity_verlet(n_p,x_vector,y_vector,z_vector,&
                               vx_vector,vy_vector,vz_vector,&
                               delta_time,mass,r_cutoff,density)

        integer(sp), intent(in)    :: n_p
        real(dp),    intent(inout) :: x_vector(n_p),y_vector(n_p),z_vector(n_p)
        real(dp),    intent(inout) :: vx_vector(n_p),vy_vector(n_p),vz_vector(n_p)
        real(dp),    intent(in)    :: delta_time,mass,r_cutoff,density
        integer(sp) :: i
        real(dp)    :: factor,force
        real(dp)    :: factor_x_old,factor_y_old,factor_z_old
        real(dp)    :: factor_x_new,factor_y_new,factor_z_new

        factor=delta_time*0.5_dp*(1._dp/mass)
        do i=1,n_p
            ! FUERZA EN EL TIEMPO ACTUAL
            force=f_lj_total(x_vector,y_vector,z_vector,r_cutoff,n_p,density)
            factor_x_old=x_vector(i)*force
            factor_y_old=y_vector(i)*force
            factor_z_old=z_vector(i)*force
            ! POSICIONES EN EL TIEMPO EVOLUCIONADO
            x_vector(i)=x_vector(i)+(vx_vector(i)+factor_x_old*factor)*delta_time
            y_vector(i)=y_vector(i)+(vy_vector(i)+factor_y_old*factor)*delta_time
            z_vector(i)=z_vector(i)+(vz_vector(i)+factor_z_old*factor)*delta_time
            call position_correction(n_p,density,x_vector(i),y_vector(i),z_vector(i))
            ! FUERZAS EN EL TIEMPO EVOLUCIONADO
            force=f_lj_total(x_vector,y_vector,z_vector,r_cutoff,n_p,density)
            factor_x_new=x_vector(i)*force
            factor_y_new=y_vector(i)*force
            factor_z_new=z_vector(i)*force
            ! VELOCIDAD EN EL TIEMPO EVOLUCIONADO
            vx_vector(i)=vx_vector(i)+(factor_x_new+factor_x_old)*factor
            vy_vector(i)=vy_vector(i)+(factor_y_new+factor_y_old)*factor
            vz_vector(i)=vz_vector(i)+(factor_z_new+factor_z_old)*factor
        end do
    end subroutine velocity_verlet

    subroutine md_initial_parameters(n_p,x_vector,y_vector,z_vector,&
                                     vx_vector,vy_vector,vz_vector,&
                                     T_adim_ref,delta_time,density)

        integer(sp), intent(in)    :: n_p ! number of particles
        real(dp),    intent(inout) :: x_vector(n_p),y_vector(n_p),z_vector(n_p)    ! position vectors
        real(dp),    intent(inout) :: vx_vector(n_p),vy_vector(n_p),vz_vector(n_p) ! velocities vectors
        real(dp),    intent(in)    :: T_adim_ref,delta_time,density ! temperature,time step,densidad
        real(dp)    :: nrand
        integer(sp) :: seed,seed_val(8),i
        real(dp)    :: kinetic_ergy_x,kinetic_ergy_y,kinetic_ergy_z ! kinetic energy
        real(dp)    :: vx_mc,vy_mc,vz_mc ! velociti center of mass
        real(dp)    :: fs_x,fs_y,fs_z    ! scale factor of the velocities

        call date_and_time(values=seed_val)
        seed=seed_val(8)*seed_val(7)*seed_val(6)+seed_val(5);call sgrnd(seed)
        do i=1,n_p ! give random velocities
            nrand=real(grnd(),dp);vx_vector(i)=(nrand-0.5_dp)
            nrand=real(grnd(),dp);vy_vector(i)=(nrand-0.5_dp)
            nrand=real(grnd(),dp);vz_vector(i)=(nrand-0.5_dp)
        end do
        ! calculamos velocidad del centro de masa
        vx_mc=sum(vx_vector(:))*(1._dp/real(n_p,dp))
        vy_mc=sum(vy_vector(:))*(1._dp/real(n_p,dp))
        vz_mc=sum(vz_vector(:))*(1._dp/real(n_p,dp))
        ! calculamos energía cinética
        kinetic_ergy_x=0.5_dp*sum(vx_vector(:)*vx_vector(:))*(1._dp/real(n_p,dp))
        kinetic_ergy_y=0.5_dp*sum(vy_vector(:)*vy_vector(:))*(1._dp/real(n_p,dp))
        kinetic_ergy_z=0.5_dp*sum(vz_vector(:)*vz_vector(:))*(1._dp/real(n_p,dp))
        ! re-escaleo de las velocidades
        fs_x=sqrt(3._dp*T_adim_ref*(1._dp/kinetic_ergy_x))
        fs_y=sqrt(3._dp*T_adim_ref*(1._dp/kinetic_ergy_y))
        fs_z=sqrt(3._dp*T_adim_ref*(1._dp/kinetic_ergy_z))
        do i=1,n_p
            ! velocity center of mass to zero
            vx_vector(i)=(vx_vector(i)-vx_mc)*fs_x
            vy_vector(i)=(vy_vector(i)-vy_mc)*fs_y
            vz_vector(i)=(vz_vector(i)-vz_mc)*fs_z
            ! position previous time step
            x_vector(i)=x_vector(i)-vx_vector(i)*delta_time
            y_vector(i)=y_vector(i)-vy_vector(i)*delta_time
            z_vector(i)=z_vector(i)-vz_vector(i)*delta_time
            ! corregimos posiciones según PBC
            call position_correction(n_p,density,x_vector(i),y_vector(i),z_vector(i))
        end do
    end subroutine md_initial_parameters
end module module_md_lennard_jones