! module of molecular dynamic to lennard jones potential
! gfortran -c module_precision.f90 module_mt19937.f90 module_md_lennard_jones.f90
module module_md_lennard_jones
    use module_precision;use module_mt19937, only: sgrnd,grnd
    implicit none
    contains

    ! FUNCIONES
    function pressure(n_p,density,mass,r_cutoff,x_vector,y_vector,z_vector,&
        vx_vector,vy_vector,vz_vector)
        integer(sp), intent(in) :: n_p
        real(dp),    intent(in) :: r_cutoff,density,mass
        real(dp),    intent(in) :: x_vector(n_p),y_vector(n_p),z_vector(n_p)
        real(dp),    intent(in) :: vx_vector(n_p),vy_vector(n_p),vz_vector(n_p)
        real(dp)                :: pressure,rij_pow02,result,T_adim,force_indiv
        integer(sp)             :: i,j
        real(dp),    parameter  :: rij_min=1.1225_dp ! 2^(1/6)*sigma
        !real(dp)                :: L,dx,dy,dz
        !L=(real(n_p,dp)*(1._dp/density))**(1._dp/3._dp)
        result=0._dp
        do j=2,n_p
            do i=1,j-1
                rij_pow02=rel_pos_correction(x_vector(i),y_vector(i),z_vector(i),&
                x_vector(j),y_vector(j),z_vector(j),n_p,density)
                cond3:  if (rij_pow02<=r_cutoff*r_cutoff) then
                            ! if (rij_pow02<=rij_min*rij_min) then
                            !     rij_pow02=rij_min*rij_min
                            !     !force_indiv=f_lj_individual(rij_pow02)
                            !     force_indiv=0._dp
                            !     exit cond3
                            ! else
                                force_indiv=f_lj_individual(rij_pow02)
                            !end if
                else
                    force_indiv=0._dp
                end if cond3
                ! dx=pbc_correction((x_vector(i)-x_vector(j)),n_p,density)
                ! dy=pbc_correction((y_vector(i)-y_vector(j)),n_p,density)
                ! dz=pbc_correction((z_vector(i)-z_vector(j)),n_p,density)
                !result=result+force_indiv*(dx*dx+dy*dy+dz*dz)
                result=result+force_indiv*rij_pow02
            end do
        end do
        T_adim=temperature(n_p,mass,vx_vector,vy_vector,vz_vector)
        pressure=density*(T_adim+(1._dp/(3._dp*real(n_p,dp)))*result)
    end function pressure

    ! compute individual lennard jones potential (simple truncation)
    function u_lj_individual(r12_pow02)
        real(dp), intent(in) :: r12_pow02 ! distancia adimensional entre pares de particulas
        real(dp)             :: r12_pow06 ! potencias de la distancia relativa
        real(dp)             :: u_lj_individual     ! adimensional lennard jones potential
        integer(sp)          :: i
        ! calculamos distancia relativa corregida según PBC
        r12_pow06=1._dp
        do i=1,3;r12_pow06=r12_pow06*r12_pow02;end do  ! (r12)^6
        u_lj_individual=4._dp*(1._dp/r12_pow06)*((1._dp/r12_pow06)-1._dp)
    end function u_lj_individual
    ! compute total lennard jones potential
    function u_lj_total(n_p,x_vector,y_vector,z_vector,r_cutoff,density)
        integer(sp), intent(in) :: n_p
        real(dp),    intent(in) :: x_vector(n_p),y_vector(n_p),z_vector(n_p)
        real(dp),    intent(in) :: r_cutoff,density
        real(dp)                :: u_lj_total,u_indiv,rij_pow02
        integer(sp)             :: i,j
        real(dp),    parameter  :: rij_min=1._dp ! sigma
        u_lj_total=0._dp
        do j=2,n_p
            do i=1,j-1
                ! calculamos distancia relativa corregida según PBC
                rij_pow02=rel_pos_correction(x_vector(i),y_vector(i),z_vector(i),&
                x_vector(j),y_vector(j),z_vector(j),n_p,density)
                cond2:  if (rij_pow02<=r_cutoff*r_cutoff) then
                            ! if (rij_pow02<=rij_min*rij_min) then
                            !     rij_pow02=rij_min*rij_min
                            !     !u_indiv=u_lj_individual(rij_pow02)
                            !     u_indiv=0._dp
                            !     exit cond2
                            ! else
                                u_indiv=u_lj_individual(rij_pow02)
                            !end if
                        else
                            u_indiv=0._dp
                end if cond2
                u_lj_total=u_lj_total+u_indiv
            end do
        end do
    end function u_lj_total

    function kinetic_ergy_total(n_p,vx_vector,vy_vector,vz_vector,mass)
        integer(sp), intent(in) :: n_p
        real(dp),    intent(in) :: vx_vector(n_p),vy_vector(n_p),vz_vector(n_p)
        real(dp),    intent(in) :: mass
        real(dp)                :: kinetic_ergy_total
        !integer(sp) :: i
        kinetic_ergy_total=0.5_dp*mass*(sum(vx_vector(:)*vx_vector(:))+&
        sum(vy_vector(:)*vy_vector(:))+sum(vz_vector(:)*vz_vector(:)))
    end function kinetic_ergy_total

    ! caclulo de la fuerza individual (par de partículas)
    function f_lj_individual(r12_pow02)
        real(dp), intent(in) :: r12_pow02   ! distancia adimensional entre pares de particulas
        real(dp)             :: r12_pow06   ! factores potencia
        real(dp)             :: f_lj_individual   ! adimensional individual lennard jones force
        integer(sp)          :: i
        r12_pow06=1._dp
        do i=1,3;r12_pow06=r12_pow06*r12_pow02;end do ! (r12)^6
        f_lj_individual=24._dp*(1._dp/(r12_pow02*r12_pow06))*(2._dp*(1._dp/r12_pow06)-1._dp)
        !f_lj_individual=-48._dp*(1._dp/r12_pow02)*((1._dp/r12_pow06)**2-0.5_dp*(1._dp/r12_pow06))
    end function f_lj_individual

    ! calculo de la componente xi de la fuerza total
    subroutine f_lj_total(x_vector,y_vector,z_vector,r_cutoff,n_p,density,&
        fx_lj_total_vector,fy_lj_total_vector,fz_lj_total_vector)
        integer(sp), intent(in)    :: n_p
        real(dp),    intent(in)    :: x_vector(n_p),y_vector(n_p),z_vector(n_p)
        real(dp),    intent(in)    :: r_cutoff,density
        real(dp),    intent(inout) :: fx_lj_total_vector(n_p),fy_lj_total_vector(n_p),&
                                      fz_lj_total_vector(n_p) ! vector de fuerzas netas
        real(dp)                   :: rij_pow02
        real(dp)                   :: force_indiv ! fuerza neta acuando en una determinada partícula
        integer(sp)                :: i,j
        real(dp),    parameter     :: rij_min=1.1225_dp ! 2^(1/6)*sigma
        real(dp)                   :: dx,dy,dz,L
        fx_lj_total_vector(:)=0._dp;fy_lj_total_vector(:)=0._dp;fz_lj_total_vector(:)=0._dp
        L=(real(n_p,dp)*(1._dp/density))**(1._dp/3._dp)
        do j=2,n_p
            do i=1,j-1
                ! calculamos distancia relativa corregida según PBC
                rij_pow02=rel_pos_correction(x_vector(i),y_vector(i),z_vector(i),&
                x_vector(j),y_vector(j),z_vector(j),n_p,density)
                cond1:  if (rij_pow02<=r_cutoff*r_cutoff) then
                            ! if (rij_pow02<=rij_min*rij_min) then
                            !     rij_pow02=rij_min*rij_min
                            !     !force_indiv=f_lj_individual(rij_pow02)
                            !     force_indiv=0._dp
                            !     exit cond1
                            ! else
                                force_indiv=f_lj_individual(rij_pow02)
                            !end if
                        else
                            force_indiv=0._dp
                end if cond1
                dx=pbc_correction((x_vector(i)-x_vector(j)),n_p,density)
                dy=pbc_correction((y_vector(i)-y_vector(j)),n_p,density)
                dz=pbc_correction((z_vector(i)-z_vector(j)),n_p,density)
                ! COMPONENTES DE LA FUERZA
                fx_lj_total_vector(i)=fx_lj_total_vector(i)+force_indiv*dx
                fx_lj_total_vector(j)=fx_lj_total_vector(j)-force_indiv*dx
                fy_lj_total_vector(i)=fy_lj_total_vector(i)+force_indiv*dy
                fy_lj_total_vector(j)=fy_lj_total_vector(j)-force_indiv*dy
                fz_lj_total_vector(i)=fz_lj_total_vector(i)+force_indiv*dz
                fz_lj_total_vector(j)=fz_lj_total_vector(j)-force_indiv*dz
            end do
        end do
    end subroutine f_lj_total

    function temperature(n_p,mass,vx_vector,vy_vector,vz_vector)
        integer(sp), intent(in) :: n_p
        real(dp),    intent(in) :: vx_vector(n_p),vy_vector(n_p),vz_vector(n_p)
        real(dp),    intent(in) :: mass
        real(dp)                :: temperature
        temperature=mass*(sum(vx_vector(:)*vx_vector(:))+&
        sum(vy_vector(:)*vy_vector(:))+sum(vz_vector(:)*vz_vector(:)))*&
        (1._dp/(3._dp*real(n_p,dp)))
    end function temperature
    
    ! corrección de las posiciones relativas (PBC)
    function rel_pos_correction(x1,y1,z1,x2,y2,z2,n_p,density)
        integer(sp), intent(in) :: n_p                ! numero total de partículas
        real(dp),    intent(in) :: density            ! densidad de partículas
        real(dp),    intent(in) :: x1,y1,z1,x2,y2,z2  ! coordenadas del par de partículas
        real(dp)                :: rel_pos_correction ! posición relativa corregida (PBC) al cuadrado
        real(dp)                :: dx,dy,dz           ! diferencia de posiciones segun coordenadas
        dx=pbc_correction((x1-x2),n_p,density)
        dy=pbc_correction((y1-y2),n_p,density)
        dz=pbc_correction((z1-z2),n_p,density)
        rel_pos_correction=dx*dx+dy*dy+dz*dz
    end function rel_pos_correction

    function pbc_correction(x,n_p,density)
        integer(sp), intent(in) :: n_p      ! numero total de partículas
        real(dp),    intent(in) :: density  ! densidad de partículas
        real(dp),    intent(in) :: x        ! variable a corregir
        real(dp)                :: pbc_correction,L
        L=(real(n_p,dp)*(1._dp/density))**(1._dp/3._dp)
        pbc_correction=(x-L*anint(x*(1._dp/L),dp)) ! MÉTODO 1
        !pbc_correction=x-L*real(int((x+0.5_dp*L)*(1._dp/L),sp),dp) ! MÉTODO 2
    end function pbc_correction

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
        x=pbc_correction(x,n_p,density)
        y=pbc_correction(y,n_p,density)
        z=pbc_correction(z,n_p,density)
    end subroutine position_correction

    ! Subroutine to set up sc and fcc lattice
    subroutine initial_lattice_configuration(n_p,density,x_vector,y_vector,z_vector,type_structure)
        integer(sp), intent(in)    :: n_p                ! numero total de partículas
        integer(sp), intent(in)    :: type_structure     ! 1 by sc lattice/ 2 by fcc lattice
        real(dp),    intent(in)    :: density            ! densidad de partículas
        real(dp),    intent(inout) :: x_vector(n_p),&    ! coordenadas de los vectores posición
                                      y_vector(n_p),&
                                      z_vector(n_p)
        integer(sp)                :: i,j,k,index2,index ! loop index
        real(dp)                   :: L                  ! logitud macroscópica por dimensión
        real(dp)                   :: a                  ! parámetro de red
        real(dp)                   :: n_unitcells        ! numero de celdas unidad por dimension
        real(dp)                   :: points_unitcells   ! número de partículas por celda unidad
        real(dp),   allocatable    :: aux_matrix(:,:)    ! matriz auxiliar con vectores primitivos (FCC)
        L=(real(n_p,dp)*(1._dp/density))**(1._dp/3._dp)
        select case (type_structure)
            case(1) ! sc laticce
                points_unitcells=1._dp
                n_unitcells=anint((real(n_p,dp)*(1._dp/points_unitcells))**(1._dp/3._dp),dp)
                a=L*(1._dp/(n_unitcells-1._dp))!a=(4*(1._dp/density))**(1._dp/3._dp)
                ! cargamos vectores de coordenadas
                index=0
                do i=1,int(n_unitcells,sp)
                    do j=1,int(n_unitcells,sp)
                        do k=1,int(n_unitcells,sp)
                            index=index+1
                            ! CENTRAMOS LA CELDA EN EL RANGO [-L/2:L/2]
                            x_vector(index)=real(i-1,dp)*a-0.5_dp*L
                            y_vector(index)=real(j-1,dp)*a-0.5_dp*L
                            z_vector(index)=real(k-1,dp)*a-0.5_dp*L
                        end do
                    end do
                end do
            case(2) ! fcc laticce
                points_unitcells=4._dp
                n_unitcells=anint((real(n_p,dp)*(1._dp/points_unitcells))**(1._dp/3._dp),dp)
                a=L*(1._dp/(n_unitcells))!a=(4*(1._dp/density))**(1._dp/3._dp)
                ! cargamos matriz con vectores primitivos (specific for FCC structure)
                allocate(aux_matrix(4,3));aux_matrix(:,:)=a*0.5_dp
                aux_matrix(1,:)=0.0_dp;aux_matrix(2,3)=0.0_dp
                aux_matrix(3,2)=0.0_dp;aux_matrix(4,1)=0.0_dp
                ! cargamos vectores de coordenadas
                index=0
                do i=1,int(n_unitcells,sp)
                    do j=1,int(n_unitcells,sp)
                        do k=1,int(n_unitcells,sp)
                            do index2=1,4
                                index=index+1
                                ! CENTRAMOS LA CELDA EN EL RANGO [-L/2:L/2]
                                x_vector(index)=(aux_matrix(index2,1)+real(i-1,dp)*a)-0.5_dp*L
                                y_vector(index)=(aux_matrix(index2,2)+real(j-1,dp)*a)-0.5_dp*L
                                z_vector(index)=(aux_matrix(index2,3)+real(k-1,dp)*a)-0.5_dp*L
                            end do
                        end do
                    end do
                end do
                deallocate(aux_matrix)
            case(3) ! fcc laticce (other method)
                points_unitcells=4._dp     
                n_unitcells=anint((real(n_p,dp)*(1._dp/points_unitcells))**(1._dp/3._dp),dp)
                a=(points_unitcells*(1._dp/density))**(1._dp/3._dp)
                index=0
                do i=0,int(n_unitcells,sp)-1
                    do j=0,int(n_unitcells,sp)-1
                        do k=0,int(n_unitcells,sp)-1
                            index=index+1
                            x_vector(index)=a*real(i,dp)-L*0.5_dp
                            y_vector(index)=a*real(j,dp)-L*0.5_dp
                            z_vector(index)=a*real(k,dp)-L*0.5_dp
                            index=index+1
                            x_vector(index)=a*(real(i,dp)+0.5_dp)-L*0.5_dp
                            y_vector(index)=a*(real(j,dp)+0.5_dp)-L*0.5_dp
                            z_vector(index)=a*real(k,dp)-L*0.5_dp
                            index=index+1
                            x_vector(index)=a*real(i,dp)-L*0.5_dp
                            y_vector(index)=a*(real(j,dp)+0.5_dp)-L*0.5_dp
                            z_vector(index)=a*(real(k,dp)+0.5_dp)-L*0.5_dp
                            index=index+1
                            x_vector(index)=a*(real(i,dp)+0.5_dp)-L*0.5_dp
                            y_vector(index)=a*real(j,dp)-L*0.5_dp
                            z_vector(index)=a*(real(k,dp)+0.5_dp)-L*0.5_dp
                        enddo
                    enddo
                enddo
            end select
    end subroutine initial_lattice_configuration

    ! SUBRUTINA DE INTEGRACIÓN DE ECUACIONES DE MOVIMIENTO
    ! subroutine velocity_verlet_v2(n_p,x_vector,y_vector,z_vector,&
    !                            vx_vector,vy_vector,vz_vector,&
    !                            delta_time,mass,r_cutoff,density,force)

    !     integer(sp), intent(in)    :: n_p
    !     real(dp),    intent(in)    :: delta_time,mass,r_cutoff,density
    !     real(dp),    intent(inout) :: x_vector(n_p),y_vector(n_p),z_vector(n_p)
    !     real(dp),    intent(inout) :: vx_vector(n_p),vy_vector(n_p),vz_vector(n_p)
    !     real(dp),    intent(inout) :: force(n_p)
    !     integer(sp)           :: i
    !     real(dp)              :: factor
    !     real(dp), allocatable :: force_x_old(:),force_y_old(:),force_z_old(:)
    !     real(dp), allocatable :: force_x_new(:),force_y_new(:),force_z_new(:)

    !     allocate(force_x_old(n_p),force_y_old(n_p),force_z_old(n_p))
    !     allocate(force_x_new(n_p),force_y_new(n_p),force_z_new(n_p))
    !     force_x_old(:)=0._dp;force_y_old(:)=0._dp;force_z_old(:)=0._dp
    !     force_x_new(:)=0._dp;force_y_new(:)=0._dp;force_z_new(:)=0._dp

    !     factor=delta_time*0.5_dp*(1._dp/mass)
    !     do i=1,n_p
    !         ! FUERZAS EN EL TIEMPO ACTUAL
    !         force_x_old(i)=x_vector(i)*force(i)
    !         force_y_old(i)=y_vector(i)*force(i)
    !         force_z_old(i)=z_vector(i)*force(i)
    !         ! POSICIONES EN EL TIEMPO EVOLUCIONADO
    !         x_vector(i)=x_vector(i)+(vx_vector(i)+force_x_old(i)*factor)*delta_time
    !         y_vector(i)=y_vector(i)+(vy_vector(i)+force_y_old(i)*factor)*delta_time
    !         z_vector(i)=z_vector(i)+(vz_vector(i)+force_z_old(i)*factor)*delta_time
    !         call position_correction(n_p,density,x_vector(i),y_vector(i),z_vector(i))
    !     end do
    !     ! actualizamos vector de fuerzas force(:)
    !     call f_lj_total(x_vector,y_vector,z_vector,r_cutoff,n_p,density,force)
    !     do i=1,n_p
    !         ! FUERZAS EN EL TIEMPO EVOLUCIONADO
    !         force_x_new(i)=x_vector(i)*force(i)
    !         force_y_new(i)=y_vector(i)*force(i)
    !         force_z_new(i)=z_vector(i)*force(i)
    !         ! VELOCIDAD EN EL TIEMPO EVOLUCIONADO
    !         vx_vector(i)=vx_vector(i)+(force_x_new(i)+force_x_old(i))*factor
    !         vy_vector(i)=vy_vector(i)+(force_y_new(i)+force_y_old(i))*factor
    !         vz_vector(i)=vz_vector(i)+(force_z_new(i)+force_z_old(i))*factor
    !     end do
    !     deallocate(force_x_old,force_y_old,force_z_old)
    !     deallocate(force_x_new,force_y_new,force_z_new)
    ! end subroutine velocity_verlet_v2

    subroutine velocity_verlet(n_p,x_vector,y_vector,z_vector,&
        vx_vector,vy_vector,vz_vector,&
        delta_time,mass,r_cutoff,density,force_x,force_y,force_z)
        integer(sp), intent(in)    :: n_p
        real(dp),    intent(inout) :: x_vector(n_p),y_vector(n_p),z_vector(n_p)
        real(dp),    intent(inout) :: vx_vector(n_p),vy_vector(n_p),vz_vector(n_p)
        real(dp),    intent(in)    :: delta_time,mass,r_cutoff,density
        real(dp),    intent(inout) :: force_x(n_p),force_y(n_p),force_z(n_p)
        integer(sp)                :: i
        real(dp)                   :: factor
        factor=delta_time*0.5_dp*(1._dp/mass)
        do i=1,n_p
            ! COMPONENTES DE LA VELOCIDAD A TIEMPO SEMI-EVOLUCIONADO
            vx_vector(i)=vx_vector(i)+force_x(i)*factor
            vy_vector(i)=vy_vector(i)+force_y(i)*factor
            vz_vector(i)=vz_vector(i)+force_z(i)*factor
            ! COMPONENTES DE LA POSICION A TIEMPO EVOLUCIONADO
            x_vector(i)=x_vector(i)+vx_vector(i)*delta_time
            y_vector(i)=y_vector(i)+vy_vector(i)*delta_time
            z_vector(i)=z_vector(i)+vz_vector(i)*delta_time
            call position_correction(n_p,density,x_vector(i),y_vector(i),z_vector(i))
        end do
        ! ACTUALIZAMOS COMPONENTES DE LA FUERZA A TIEMPO EVOLUCIONADO
        call f_lj_total(x_vector,y_vector,z_vector,r_cutoff,n_p,density,force_x,force_y,force_z)
        do i=1,n_p
            ! COMPONENTES DE LA VELOCIDAD A TIEMPO EVOLUCIONADO
            vx_vector(i)=vx_vector(i)+force_x(i)*factor
            vy_vector(i)=vy_vector(i)+force_y(i)*factor
            vz_vector(i)=vz_vector(i)+force_z(i)*factor
        end do
    end subroutine velocity_verlet

    subroutine md_initial_parameters(n_p,x_vector,y_vector,z_vector,&
                                     vx_vector,vy_vector,vz_vector,&
                                     T_adim_ref,delta_time,density,mass)

        integer(sp), intent(in)    :: n_p ! number of particles
        real(dp),    intent(inout) :: x_vector(n_p),y_vector(n_p),z_vector(n_p)    ! position vectors
        real(dp),    intent(inout) :: vx_vector(n_p),vy_vector(n_p),vz_vector(n_p) ! velocities vectors
        real(dp),    intent(in)    :: T_adim_ref,delta_time,density,mass ! temperature,time step,densidad,masa
        real(dp)    :: nrand
        integer(sp) :: seed,seed_val(8),i
        real(dp)    :: vx_mc,vy_mc,vz_mc ! velocity center of mass

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
        ! velocity center of mass to zero and rescaling
        vx_vector(:)=vx_vector(:)-vx_mc
        vy_vector(:)=vy_vector(:)-vy_mc
        vz_vector(:)=vz_vector(:)-vz_mc
        call rescaling_velocities(n_p,vx_vector,vy_vector,vz_vector,T_adim_ref,mass)
        ! do i=1,n_p
        !     ! position previous time step
        !     x_vector(i)=x_vector(i)-vx_vector(i)*delta_time
        !     y_vector(i)=y_vector(i)-vy_vector(i)*delta_time
        !     z_vector(i)=z_vector(i)-vz_vector(i)*delta_time
        !     ! corregimos posiciones según PBC
        !     call position_correction(n_p,density,x_vector(i),y_vector(i),z_vector(i))
        ! end do
    end subroutine md_initial_parameters
end module module_md_lennard_jones