! make clean && make lennard_jones_fluid.o && ./lennard_jones_fluid.o
program lennard_jones_fluid
    use module_precision
    implicit none

    real(dp)                 :: density          ! densidad de partículas
    integer(sp)              :: n_p              ! numero total de partículas (125,216,512)
    integer(sp)              :: i,istat
    real(dp),    allocatable :: sc_vector_x(:),&
                                sc_vector_y(:),&
                                sc_vector_z(:)

    20 format(2(E12.4,x),E12.4);21 format(2(A12,x),A12)

    ! estructura cubica simple (SC)
    n_p=125_sp;density=1._dp
    allocate(sc_vector_x(n_p),sc_vector_y(n_p),sc_vector_z(n_p))
    sc_vector_x=0._dp;sc_vector_y=0._dp;sc_vector_z=0._dp
    call initial_particle_configuration_sc(n_p,density,sc_vector_x,sc_vector_y,sc_vector_z)
    open(10,file='../results/lennard_jones_fluid_01.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(10file) = ',istat
    write(10,21) 'rx_sc','ry_sc','rz_sc'
    do i=1,n_p;write(10,20) sc_vector_x(i),sc_vector_y(i),sc_vector_z(i);end do
    close(10)
    deallocate(sc_vector_x,sc_vector_y,sc_vector_z)

    ! estructura cubica centrada en el cuerpo (FCC)
    n_p=256_sp;density=0.8_dp
    allocate(sc_vector_x(n_p),sc_vector_y(n_p),sc_vector_z(n_p))
    sc_vector_x=0._dp;sc_vector_y=0._dp;sc_vector_z=0._dp
    call initial_particle_configuration_fcc(n_p,density,sc_vector_x,sc_vector_y,sc_vector_z)
    open(10,file='../results/lennard_jones_fluid_02.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(10file) = ',istat
    write(10,21) 'rx_fcc','ry_fcc','rz_fcc'
    do i=1,n_p;write(10,20) sc_vector_x(i),sc_vector_y(i),sc_vector_z(i);end do
    close(10)
    deallocate(sc_vector_x,sc_vector_y,sc_vector_z)

end program lennard_jones_fluid

subroutine initial_particle_configuration_sc(n_p,density,x_coord_vector,y_coord_vector,z_coord_vector)
    use module_precision
    implicit none
    integer(sp), intent(in)    :: n_p                   ! numero total de partículas
    real(dp),    intent(in)    :: density               ! densidad de partículas
    real(dp),    intent(inout) :: x_coord_vector(n_p),& ! coordenadas de los vectores posición
                                  y_coord_vector(n_p),&
                                  z_coord_vector(n_p)

    integer(sp)           :: i,j,k,index
    integer(sp)           :: m                  ! numero total de coordenadas (x,y,z)
    real(dp)              :: L                  ! longitud macroscópica
    real(dp)              :: a                  ! parámetro de red
    real(dp), allocatable :: aux_vector_x(:),&
                             aux_vector_y(:),&
                             aux_vector_z(:)

    L=(real(n_p,dp)*(1._dp/density))**(1._dp/3._dp)

    m=int(anint(n_p**(1._dp/3._dp),dp),sp)
    a=L*(1._dp/real(m-1,dp))
    allocate(aux_vector_x(m),aux_vector_y(m),aux_vector_z(m))
    aux_vector_x=0._dp;aux_vector_y=0._dp;aux_vector_z=0._dp
    ! creamos una configuración SC
    ! llenamos los vectores por coordenadas (todos diferentes)
    do i=1,m
        aux_vector_x(i)=real(i-1,dp)*a
        aux_vector_y(i)=real(i-1,dp)*a
        aux_vector_z(i)=real(i-1,dp)*a
    end do
    ! llenamos las coordenadas de cada particula
    index=0
    do i=1,m;do j=1,m;do k=1,m
        index=index+1
        x_coord_vector(index)=aux_vector_x(i)
        y_coord_vector(index)=aux_vector_y(j)
        z_coord_vector(index)=aux_vector_z(k)
    end do;end do;end do
    deallocate(aux_vector_x,aux_vector_y,aux_vector_z)
end subroutine initial_particle_configuration_sc

! Subroutine to set up fcc lattice
subroutine initial_particle_configuration_fcc(n_p,density,x_coord_vector,y_coord_vector,z_coord_vector)
    use module_precision
    implicit none
    integer(sp), intent(in)    :: n_p                   ! numero total de partículas
    real(dp),    intent(in)    :: density               ! densidad de partículas
    real(dp),    intent(inout) :: x_coord_vector(n_p),& ! coordenadas de los vectores posición
                                  y_coord_vector(n_p),&
                                  z_coord_vector(n_p)

    integer(sp)           :: i,j,k,index2,index
    real(dp)              :: L                  ! logitud macroscópica por dimensión
    real(dp)              :: a                  ! parámetro de red
    real(dp)              :: n_unitcells ! numero de celdas unidad por dimension
    real(dp), allocatable :: aux_matrix(:,:)

    allocate(aux_matrix(4,3));aux_matrix(4,3)=0._dp
    
    L=(n_p*(1._dp/density))**(1._dp/3._dp)   
    n_unitcells=anint((n_p*0.25_dp)**(1._dp/3._dp),dp) 
    a=L*(1._dp/n_unitcells)
    !a=(4*(1._dp/density))**(1._dp/3._dp)

    aux_matrix(:,:)=a*0.5_dp
    aux_matrix(1,:)=0.0_dp;aux_matrix(2,3)=0.0_dp
    aux_matrix(3,2)=0.0_dp;aux_matrix(4,1)=0.0_dp
    
    ! Going through all boxes in the system
    index=0
    do i=1,int(n_unitcells,sp);do j=1,int(n_unitcells,sp);do k=1,int(n_unitcells,sp);do index2=1,4
        index=index+1
        x_coord_vector(index) = aux_matrix(index2,1) + real(i-1,dp)*a
        y_coord_vector(index) = aux_matrix(index2,2) + real(j-1,dp)*a
        z_coord_vector(index) = aux_matrix(index2,3) + real(k-1,dp)*a
    end do;end do;end do;end do

    deallocate(aux_matrix)
    
end subroutine initial_particle_configuration_fcc

function u_lj(x,y,z)
    use module_precision
    implicit none
    real(dp), intent(in) :: x,y,z           ! coordenadas
    real(dp)             :: factor1,factor2 ! factores posiciones
    real(dp)             :: r               ! vector posición adimensional
    real(dp)             :: u_lj            ! lennard_jones_potential
    integer(sp)          :: i

    r=sqrt(x*x+y*y+z*z)
    factor1=1._dp;factor2=1._dp
    do i=1,12;factor1=factor1*(1._dp/r);end do
    do i=1,6;factor2=factor2*(1._dp/r);end do

        u_lj=4*(factor1-factor2)
end function u_lj