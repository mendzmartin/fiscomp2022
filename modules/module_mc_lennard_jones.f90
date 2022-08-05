! Module of monte carlo dynamic (metropolis method) to lennard jones potential
! Compile command
!   gfortran -c module_precision.f90 module_mt19937.f90 module_mc_lennard_jones.f90
module module_mc_lennard_jones
    use module_precision;use module_mt19937, only: sgrnd,grnd
    implicit none
    contains

    ! FUNCIONES
    ! FUNCION PARA CALCULAR EL FACTOR DE ESTRUCTURA ESTÁTICO (PARAMETRO DE ORDEN CRISTALINO)
    function static_structure_factor(n_p,density,x_vector,y_vector,z_vector)
        integer(sp), intent(in) :: n_p                                       ! numero de partículas
        real(dp),    intent(in) :: density                                   ! densidad de partículas
        real(dp),    intent(in) :: x_vector(n_p),y_vector(n_p),z_vector(n_p) ! coordenadas del vector posicion
        real(dp)                :: static_structure_factor                   ! funcion de estructura
        real(dp),    parameter  :: pi=4._dp*atan(1._dp)
        real(dp)                :: a                                         ! parámetro de red
        real(dp)                :: points_unitcells                          ! número de partículas por celda unidad
        real(dp)                :: factor1,factor2
        real(dp)                :: kx,ky,kz                                  ! componentes del vector de onda
        integer(sp)             :: i
    
        ! valido únicamente para una red FCC
        points_unitcells=4._dp
        a=(points_unitcells*(1._dp/density))**(1._dp/3._dp)
        factor1=0._dp;factor2=0._dp
        do i=1,n_p
            ! cambiar en caso de tener otro vector de onda
            kx=-2._dp*pi*(1._dp/a)*x_vector(i)
            ky=2._dp*pi*(1._dp/a)*y_vector(i)
            kz=-2._dp*pi*(1._dp/a)*z_vector(i)
            factor1=factor1+cos(kx+ky+kz)
            factor2=factor2+sin(kx+ky+kz)
        end do
        factor1=factor1*factor1;factor2=factor2*factor2
        static_structure_factor=(1._dp/real(n_p*n_p,dp))*(factor1+factor2)
    end function static_structure_factor

    function osmotic_pressure(n_p,density,mass,r_cutoff,x_vector,y_vector,z_vector)
        integer(sp), intent(in) :: n_p
        real(dp),    intent(in) :: r_cutoff,density,mass
        real(dp),    intent(in) :: x_vector(n_p),y_vector(n_p),z_vector(n_p)
        real(dp)                :: osmotic_pressure,rij_pow02,result,force_indiv
        integer(sp)             :: i,j
        result=0._dp
        do j=2,n_p
            do i=1,j-1
                rij_pow02=rel_pos_correction(x_vector(i),y_vector(i),z_vector(i),&
                x_vector(j),y_vector(j),z_vector(j),n_p,density)
                if (rij_pow02<=r_cutoff*r_cutoff) then
                    force_indiv=f_lj_individual(rij_pow02)
                else; force_indiv=0._dp; end if
                result=result+force_indiv*rij_pow02
            end do
        end do
        osmotic_pressure=density*(1._dp/(3._dp*real(n_p,dp)))*result
    end function osmotic_pressure

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
    ! compute total lennard jones potential (TRUNCADO Y DESPLAZADO)
    function u_lj_total(n_p,x_vector,y_vector,z_vector,r_cutoff,density)
        integer(sp), intent(in) :: n_p
        real(dp),    intent(in) :: x_vector(n_p),y_vector(n_p),z_vector(n_p)
        real(dp),    intent(in) :: r_cutoff,density
        real(dp)                :: u_lj_total,u_indiv,rij_pow02
        integer(sp)             :: i,j
        u_lj_total=0._dp
        do j=2,n_p
            do i=1,j-1
                ! calculamos distancia relativa corregida según PBC
                rij_pow02=rel_pos_correction(x_vector(i),y_vector(i),z_vector(i),&
                x_vector(j),y_vector(j),z_vector(j),n_p,density)
                if (rij_pow02<=r_cutoff*r_cutoff) then
                    u_indiv=u_lj_individual(rij_pow02)-u_lj_individual(r_cutoff*r_cutoff)
                else; u_indiv=0._dp; end if
                u_lj_total=u_lj_total+u_indiv
            end do
        end do
        !u_lj_total=u_lj_total*(1._dp/real(n_p,dp))
    end function u_lj_total

    ! compute total lennard jones potential (TRUNCADO Y DESPLAZADO)
    function u_lj_reduced(n_p,x_vector,y_vector,z_vector,r_cutoff,density,index)
        integer(sp), intent(in) :: n_p,index
        real(dp),    intent(in) :: x_vector(n_p),y_vector(n_p),z_vector(n_p)
        real(dp),    intent(in) :: r_cutoff,density
        real(dp)                :: u_lj_reduced,u_indiv,rij_pow02
        integer(sp)             :: i,j ! indices para etiquetar par de partículas
        
        u_lj_reduced=0._dp
        j=index
        do i=1,index-1
            ! calculamos distancia relativa corregida según PBC
            rij_pow02=rel_pos_correction(x_vector(i),y_vector(i),z_vector(i),&
                x_vector(j),y_vector(j),z_vector(j),n_p,density)
            if (rij_pow02<=r_cutoff*r_cutoff) then
                u_indiv=u_lj_individual(rij_pow02)-u_lj_individual(r_cutoff*r_cutoff)
            else; u_indiv=0._dp; end if
            u_lj_reduced=u_lj_reduced+u_indiv
        end do

        i=index
        do j=index+1,n_p
            ! calculamos distancia relativa corregida según PBC
            rij_pow02=rel_pos_correction(x_vector(i),y_vector(i),z_vector(i),&
                x_vector(j),y_vector(j),z_vector(j),n_p,density)
            if (rij_pow02<=r_cutoff*r_cutoff) then
                u_indiv=u_lj_individual(rij_pow02)-u_lj_individual(r_cutoff*r_cutoff)
            else; u_indiv=0._dp; end if
            u_lj_reduced=u_lj_reduced+u_indiv
        end do
    end function u_lj_reduced

    function delta_u_lj(n_p,x_vector,y_vector,z_vector,r_cutoff,density,x_value,y_value,z_value,index)
        integer(sp), intent(in)    :: n_p,index ! numero de partículas e indície de la partícula desplazada
        real(dp),    intent(inout) :: x_vector(n_p),y_vector(n_p),z_vector(n_p)
        real(dp),    intent(in)    :: x_value,y_value,z_value ! valores posición originales
        real(dp),    intent(in)    :: r_cutoff,density
        real(dp)                   :: delta_u_lj
        real(dp)                   :: x_value_new,y_value_new,z_value_new

        ! guardamos las posiciones desplazadas
        x_value_new=x_vector(index);y_value_new=y_vector(index);z_value_new=z_vector(index)

        ! computamos energía con posiciones sin desplazar (valores originales)
        x_vector(index)=x_value;y_vector(index)=y_value;z_vector(index)=z_value
        delta_u_lj=-u_lj_reduced(n_p,x_vector,y_vector,z_vector,r_cutoff,density,index)
        ! computamos energía con posiciones desplazadas (valores originales)
        x_vector(index)=x_value_new;y_vector(index)=y_value_new;z_vector(index)=z_value_new
        delta_u_lj=delta_u_lj+u_lj_reduced(n_p,x_vector,y_vector,z_vector,r_cutoff,density,index)
        delta_u_lj=delta_u_lj*(0.5_dp/real(n_p,dp))
    end function delta_u_lj

    ! caclulo de la fuerza individual (par de partículas)
    function f_lj_individual(r12_pow02)
        real(dp), intent(in) :: r12_pow02   ! distancia adimensional entre pares de particulas
        real(dp)             :: r12_pow06   ! factores potencia
        real(dp)             :: f_lj_individual   ! adimensional individual lennard jones force
        integer(sp)          :: i
        r12_pow06=1._dp
        do i=1,3;r12_pow06=r12_pow06*r12_pow02;end do ! (r12)^6
        f_lj_individual=24._dp*(1._dp/(r12_pow02*r12_pow06))*(2._dp*(1._dp/r12_pow06)-1._dp)
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
        real(dp)                   :: dx,dy,dz,L
        fx_lj_total_vector(:)=0._dp;fy_lj_total_vector(:)=0._dp;fz_lj_total_vector(:)=0._dp
        L=(real(n_p,dp)*(1._dp/density))**(1._dp/3._dp)
        do j=2,n_p
            do i=1,j-1
                ! calculamos distancia relativa corregida según PBC
                rij_pow02=rel_pos_correction(x_vector(i),y_vector(i),z_vector(i),&
                x_vector(j),y_vector(j),z_vector(j),n_p,density)

                if (rij_pow02<=r_cutoff*r_cutoff) then
                    force_indiv=f_lj_individual(rij_pow02)
                else; force_indiv=0._dp;end if

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

    subroutine f_lj_total_linkedlist(x_vector,y_vector,z_vector,r_cutoff,n_p,density,&
        fx_lj_total_vector,fy_lj_total_vector,fz_lj_total_vector,m,map,list,head)
        integer(sp), intent(in)    :: n_p,m,map(13*m*m*m)
        real(dp),    intent(in)    :: x_vector(n_p),y_vector(n_p),z_vector(n_p)
        real(dp),    intent(in)    :: r_cutoff,density
        real(dp),    intent(inout) :: fx_lj_total_vector(n_p),fy_lj_total_vector(n_p),&
                                      fz_lj_total_vector(n_p) ! vector de fuerzas netas
        integer(sp), intent(inout) :: list(n_p),head(m*m*m)
        real(dp)                   :: rij_pow02
        ! fuerza neta actuando en una determinada partícula
        real(dp)                   :: force_indiv
        integer(sp)                :: i,j
        real(dp)                   :: dx,dy,dz,L
        integer(sp)                :: icell,jcell,jcell0,nabor

        fx_lj_total_vector(:)=0._dp;fy_lj_total_vector(:)=0._dp;fz_lj_total_vector(:)=0._dp
        L=(real(n_p,dp)*(1._dp/density))**(1._dp/3._dp)

        do icell=1,m*m*m ! celdas
            i=head(icell)
            do while (i/=0)
                j=list(i)
                do while (j/=0) ! en la celda
                    rij_pow02=rel_pos_correction(x_vector(i),y_vector(i),z_vector(i),&
                        x_vector(j),y_vector(j),z_vector(j),n_p,density)
                    if (rij_pow02==0._dp) stop
                    if (rij_pow02<=r_cutoff*r_cutoff) then
                        force_indiv=f_lj_individual(rij_pow02)
                    else; force_indiv=0._dp;end if
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
                    ! actualizo el indice que recorre list(:)
                    j = list(j)
                end do
                jcell0=13*(icell-1)
                do nabor=1,13 ! celdas vecinas
                    jcell=map(jcell0+nabor)
                    j=head(jcell)
                    do while (j/=0) ! en la celda
                        rij_pow02=rel_pos_correction(x_vector(i),y_vector(i),z_vector(i),&
                            x_vector(j),y_vector(j),z_vector(j),n_p,density)
                        if (rij_pow02==0._dp) stop
                        if (rij_pow02<=r_cutoff*r_cutoff) then
                            force_indiv=f_lj_individual(rij_pow02)
                        else; force_indiv=0._dp;end if
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
                        j=list(j)
                    end do
                end do
                i=list(i)
            end do
        end do ! celdas   
    end subroutine f_lj_total_linkedlist

    subroutine maps(m,map)
        integer(sp), intent(in)    :: m
        integer(sp), intent(inout) :: map(13*m*m*m)
        integer(sp)                :: index_x,index_y,index_z,imap
        do index_x=1,m;do index_y=1,m;do index_z=1,m
            imap=(icell_function(index_x,index_y,index_z,m)-1)*13
            map(imap+1)=icell_function(index_x+1,index_y,index_z,M)
            map(imap+2)=icell_function(index_x+1,index_y+1,index_z,M)
            map(imap+3)=icell_function(index_x,index_y+1,index_z,M)
            map(imap+4)=icell_function(index_x-1,index_y+1, index_z,M)
            map(imap+5)=icell_function(index_x+1,index_y,index_z-1,M)
            map(imap+6)=icell_function(index_x+1,index_y+1,index_z-1,M)
            map(imap+7)=icell_function(index_x,index_y+1,index_z-1,M)
            map(imap+8)=icell_function(index_x-1,index_y+1,index_z-1,M)
            map(imap+9)=icell_function(index_x+1,index_y,index_z+1,M)
            map(imap+10)=icell_function(index_x+1,index_y+1,index_z+1,M)
            map(imap+11)=icell_function(index_x,index_y+1,index_z+1,M)
            map(imap+12)=icell_function(index_x-1,index_y+1,index_z+1,M)
            map(imap+13)=icell_function(index_x,index_y,index_z+1,M)
        end do;end do;end do
    end subroutine maps

    function icell_function(index_x,index_y,index_z,m)
        integer(sp), intent(in) :: index_x,index_y,index_z,m
        integer(sp) :: icell_function
        icell_function=1+mod(index_x-1+m,m)+mod(index_y-1+m,m)*m+mod(index_z-1+m,m)*m*m
    end function icell_function
    
    subroutine links(n_p,m,L,head,list,x_vector,y_vector,z_vector)
        integer(sp), intent(in)    :: m,n_p
        real(dp),    intent(in)    :: L,x_vector(n_p),y_vector(n_p),z_vector(n_p)
        integer(sp), intent(inout) :: head(m*m*m),list(n_p)
        real(dp)                   :: Lc_inv  ! inversa del lado de la celda
        integer(sp)                :: i,icell 
        ! celda de -L/2 a L/2
        head(:)=0_sp;Lc_inv=real(m,dp)*(1._dp/L)
        do i=1,n_p
            icell=1+int((x_vector(i)+0.5_dp*L)*Lc_inv,sp)+&
                    int((y_vector(i)+0.5_dp*L)*Lc_inv,sp)*m+&
                    int((z_vector(i)+0.5_dp*L)*Lc_inv,sp)*m*m
            list(i)=head(icell);head(icell)=i
        end do
    end subroutine links
    
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
        pbc_correction=(x-L*anint(x*(1._dp/L),dp))
    end function pbc_correction

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
            case(1) ! sc laticce (1st method)
                points_unitcells=1._dp
                n_unitcells=anint((real(n_p,dp)*(1._dp/points_unitcells))**(1._dp/3._dp),dp)
                a=L*(1._dp/(n_unitcells-1._dp))!a=(4*(1._dp/density))**(1._dp/3._dp)
                ! cargamos vectores de coordenadas
                index=0
                do i=1,int(n_unitcells,sp);do j=1,int(n_unitcells,sp);do k=1,int(n_unitcells,sp)
                    index=index+1
                    ! CENTRAMOS LA CELDA EN EL RANGO [-L/2:L/2]
                    x_vector(index)=real(i-1,dp)*a-0.5_dp*L
                    y_vector(index)=real(j-1,dp)*a-0.5_dp*L
                    z_vector(index)=real(k-1,dp)*a-0.5_dp*L
                end do;end do;end do
            case(2) ! fcc laticce (2nd method)
                points_unitcells=4._dp
                n_unitcells=anint((real(n_p,dp)*(1._dp/points_unitcells))**(1._dp/3._dp),dp)
                a=L*(1._dp/(n_unitcells))!a=(4*(1._dp/density))**(1._dp/3._dp)
                ! cargamos matriz con vectores primitivos (specific for FCC structure)
                allocate(aux_matrix(4,3));aux_matrix(:,:)=a*0.5_dp
                aux_matrix(1,:)=0.0_dp;aux_matrix(2,3)=0.0_dp
                aux_matrix(3,2)=0.0_dp;aux_matrix(4,1)=0.0_dp
                ! cargamos vectores de coordenadas
                index=0
                do i=1,int(n_unitcells,sp);do j=1,int(n_unitcells,sp);do k=1,int(n_unitcells,sp)
                    do index2=1,4
                        index=index+1
                        ! CENTRAMOS LA CELDA EN EL RANGO [-L/2:L/2]
                        x_vector(index)=(aux_matrix(index2,1)+real(i-1,dp)*a)-0.5_dp*L
                        y_vector(index)=(aux_matrix(index2,2)+real(j-1,dp)*a)-0.5_dp*L
                        z_vector(index)=(aux_matrix(index2,3)+real(k-1,dp)*a)-0.5_dp*L
                    end do
                end do;end do;end do
                deallocate(aux_matrix)
            case(3) ! fcc laticce (3rd method)
                points_unitcells=4._dp     
                n_unitcells=anint((real(n_p,dp)*(1._dp/points_unitcells))**(1._dp/3._dp),dp)
                a=(points_unitcells*(1._dp/density))**(1._dp/3._dp)
                index=0
                do i=0,int(n_unitcells,sp)-1;do j=0,int(n_unitcells,sp)-1;do k=0,int(n_unitcells,sp)-1
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
                end do;end do;end do
            end select
    end subroutine initial_lattice_configuration

    ! SUBRUTINA PARA AJUSTE DE DESPLAZAMIENTO MÁXIMO (optimizado)
    subroutine max_displacement_adjusting(n_p,x_vector,y_vector,z_vector,&
        T_adim,r_cutoff,density,delta_x,delta_y,delta_z)
        integer(sp), intent(in)    :: n_p                                    ! cantidad de partículas
        real(dp),    intent(inout) :: x_vector(n_p),y_vector(n_p),z_vector(n_p) ! vectores posicion
        real(dp),    intent(in)    :: T_adim                                 ! Temperatura de equilibrio
        real(dp),    intent(in)    :: r_cutoff ! distancia de truncado
        real(dp),    intent(in)    :: density
        real(dp),    intent(inout) :: delta_x,delta_y,delta_z
        real(dp)    :: deltaU_adim,U_adim_old ! variación de energía
        real(dp)    :: x_old,y_old,z_old ! posiciones sin deplazar con PBC
        real(dp)    :: L                 ! desplazamientos
        real(dp)    :: nrand             ! numero random
        integer(sp) :: seed,seed_val(8),MC_index,index
        integer(sp) :: counter           ! contador de aceptaciones
        real(dp)    :: accept_prob       ! acceptance_probability
        logical     :: end_loop           
        
        L=(real(n_p,dp)*(1._dp/density))**(1._dp/3._dp)

        call date_and_time(values=seed_val)
        seed=seed_val(8)*seed_val(7)*seed_val(6)+seed_val(5);call sgrnd(seed)

        end_loop=.false.
        do while (end_loop.eqv..false.)
            counter=0_sp
            ! hago un paso de Monte Carlo (MC_step)
            do MC_index=1,n_p
                ! random integer from 1 to n_p (elijo particulas al azar y desplazo)
                nrand=real(grnd(),dp);index=1_sp+floor(n_p*nrand,sp)
                ! guardamos posicion en tiempo actual
                x_old=x_vector(index);y_old=y_vector(index);z_old=z_vector(index)
                ! evolucionamos posiciones con y sin PBC
                nrand=real(grnd(),dp);x_vector(index)=x_old+(nrand-0.5_dp)*delta_x
                nrand=real(grnd(),dp);y_vector(index)=y_old+(nrand-0.5_dp)*delta_y
                nrand=real(grnd(),dp);z_vector(index)=z_old+(nrand-0.5_dp)*delta_z
                ! corrección de posiciones según PBC
                call position_correction(n_p,density,x_vector(index),y_vector(index),z_vector(index))
                ! variación de energía interna
                deltaU_adim=delta_u_lj(n_p,x_vector,y_vector,z_vector,r_cutoff,density,x_old,y_old,z_old,index)
                cond1:  if (deltaU_adim<0._dp) then;counter=counter+1_sp;exit cond1
                    else
                        if (T_adim==0._dp) exit cond1
                        nrand=real(grnd(),dp)
                        if ((deltaU_adim*(1._dp/T_adim))<abs(log(tiny(1._dp)))) then
                            ! método de metrópolis
                            if (exp(-deltaU_adim*(1._dp/T_adim))>nrand) then;counter=counter+1_sp;end if
                            exit cond1
                        end if
                end if cond1
                x_vector(index)=x_old;y_vector(index)=y_old;z_vector(index)=z_old
            end do
            ! evaluamos según probabilidad de aceptación
            accept_prob=real(counter,dp)/real(n_p,dp)
            if(accept_prob>0.52_dp) then
                delta_x=delta_x*1.05_dp;delta_y=delta_y*1.05_dp;delta_z=delta_z*1.05_dp
            else if(accept_prob<0.48_dp) then
                delta_x=delta_x*0.95_dp;delta_y=delta_y*0.95_dp;delta_z=delta_z*0.95_dp
            else;end_loop=.true.;end if
        end do
    end subroutine max_displacement_adjusting

    ! SUBRUTINA DE INTEGRACIÓN DE ECUACIONES DE MOVIMIENTO
    subroutine evolution_monte_carlo(n_p,x_vector,y_vector,z_vector,&
        x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,&
        U_adim,T_adim,r_cutoff,density,delta_x,delta_y,delta_z)
        integer(sp), intent(in)    :: n_p                                       ! cantidad de partículas
        real(dp),    intent(inout) :: x_vector(n_p),y_vector(n_p),z_vector(n_p) ! vectores posicion
        real(dp),    intent(inout) :: x_vector_noPBC(n_p),y_vector_noPBC(n_p),&
                                      z_vector_noPBC(n_p)                       ! componentes del vector posición sin PBC
        real(dp),    intent(inout) :: U_adim                                    ! Energía interna
        real(dp),    intent(in)    :: T_adim                                    ! Temperatura de equilibrio
        real(dp),    intent(in)    :: r_cutoff ! distancia de truncado
        real(dp),    intent(in)    :: density
        real(dp),    intent(in)    :: delta_x,delta_y,delta_z                   ! desplazamientos optimizados
        real(dp)    :: deltaU_adim                          ! variación de energía
        real(dp)    :: x_old,y_old,z_old                    ! posiciones sin deplazar
        real(dp)    :: x_old_noPBC,y_old_noPBC,z_old_noPBC
        real(dp)    :: L
        real(dp)    :: nrand                                ! numero random
        integer(sp) :: seed,seed_val(8),MC_index,index
        real(dp)    :: U_adim_old,U_adim_new,deltaU_adim_good

        L=(real(n_p,dp)*(1._dp/density))**(1._dp/3._dp)

        call date_and_time(values=seed_val)
        seed=seed_val(8)*seed_val(7)*seed_val(6)+seed_val(5);call sgrnd(seed)

        ! hago un paso de Monte Carlo (MC_step)
        do MC_index=1,n_p
            ! random integer from 1 to n_p (elijo particulas al azar y desplazo)
            nrand=real(grnd(),dp);index=1_sp+floor(n_p*nrand,sp)
            ! guardamos posicion en tiempo actual
            x_old=x_vector(index);x_old_noPBC=x_vector_noPBC(index)
            y_old=y_vector(index);y_old_noPBC=y_vector_noPBC(index)
            z_old=z_vector(index);z_old_noPBC=z_vector_noPBC(index)
            ! evolucionamos posiciones con y sin PBC
            nrand=real(grnd(),dp);x_vector(index)=x_old+(nrand-0.5_dp)*delta_x
            x_vector_noPBC(index)=x_old_noPBC+(nrand-0.5_dp)*delta_x
            nrand=real(grnd(),dp);y_vector(index)=y_old+(nrand-0.5_dp)*delta_y
            y_vector_noPBC(index)=y_old_noPBC+(nrand-0.5_dp)*delta_y
            nrand=real(grnd(),dp);z_vector(index)=z_old+(nrand-0.5_dp)*delta_z
            z_vector_noPBC(index)=z_old_noPBC+(nrand-0.5_dp)*delta_z
            ! corrección de posiciones según PBC
            call position_correction(n_p,density,x_vector(index),y_vector(index),z_vector(index))
            ! variación de energía interna
            deltaU_adim=delta_u_lj(n_p,x_vector,y_vector,z_vector,r_cutoff,density,x_old,y_old,z_old,index)
            !U_adim_old=U_adim;U_adim_new=u_lj_total(n_p,x_vector,y_vector,z_vector,r_cutoff,density)
            !deltaU_adim=U_adim_new-U_adim_old
            !write(*,*) MC_index,deltaU_adim,deltaU_adim_good
            cond1:  if (deltaU_adim<0._dp) then;U_adim=U_adim+deltaU_adim;exit cond1
                else
                    if (T_adim==0._dp) exit cond1
                    nrand=real(grnd(),dp)
                    if ((deltaU_adim*(1._dp/T_adim))<abs(log(tiny(1._dp)))) then
                        ! método de metrópolis
                        if (exp(-deltaU_adim*(1._dp/T_adim))>nrand) U_adim=U_adim+deltaU_adim
                        exit cond1
                    else
                        x_vector(index)=x_old;x_vector_noPBC(index)=x_old_noPBC
                        y_vector(index)=y_old;y_vector_noPBC(index)=y_old_noPBC
                        z_vector(index)=z_old;z_vector_noPBC(index)=z_old_noPBC
                    end if
            end if cond1
            ! x_vector(index)=x_old;x_vector_noPBC(index)=x_old_noPBC
            ! y_vector(index)=y_old;y_vector_noPBC(index)=y_old_noPBC
            ! z_vector(index)=z_old;z_vector_noPBC(index)=z_old_noPBC
            ! U_adim=U_adim_old
        end do
    end subroutine evolution_monte_carlo

    ! SUBRUTINA PARA CALCULAR LA FUNCION DE DISTRIBUCIÓN RADIAL (FUNCION DE CORRELACIÓN)
    subroutine radial_ditribution_function(file_name,n_p,density,x_vector,y_vector,z_vector,n_bins,g)
        character(len=*), intent(in)    :: file_name                                 ! nombre del archivo de datos
        integer(sp),      intent(in)    :: n_p,n_bins                                ! numero de partículas y numero total de bins
        real(dp),         intent(in)    :: x_vector(n_p),y_vector(n_p),z_vector(n_p) ! coordenadas del vector posicion
        real(dp),         intent(in)    :: density                                   ! densidad de partículas
        real(dp),         intent(inout) :: g(n_bins)                                 ! funcion distribución

        integer(sp) :: ngr              ! contador de particulas
        real(dp)    :: rij_pow02,rij    ! distancia relativa entre par de partículas
        real(dp)    :: L                ! logitud macroscópica
        real(dp)    :: step_bins        ! tamaño de bins 
        real(dp)    :: volume           ! volumen de particulas
        real(dp)    :: nid              ! numero de partículas del gas ideal en volume
        integer(sp) :: i,j,index,istat  ! loop variables

        ! initialization
        L=(real(n_p,dp)*(1._dp/density))**(1._dp/3._dp)
        step_bins=L/(2*n_bins)
        g(:)=0._dp

        ! sample
        ngr=0_sp;ngr=ngr+1_sp
        do j=2,n_p;do i=1,j-1
            ! relative distances (pow 2) with pbc correction
            rij_pow02=rel_pos_correction(x_vector(i),y_vector(i),z_vector(i),&
            x_vector(j),y_vector(j),z_vector(j),n_p,density)
            ! relative distances
            rij=sqrt(rij_pow02)
            if (rij<=L/2) then ! only within half the box length
                index=int(rij/step_bins)
                g(index)=g(index)+2._dp ! contribution for particle i and j
            end if
        end do;end do
        
        ! result initialization (determine g(rij))
        open(100,file=file_name,status='replace',action='write',iostat=istat)
            if (istat/=0) write(*,*) 'ERROR! istat(11file) = ',istat
            101 format(E12.4,x,E12.4);102 format(A12,x,A12)
            write(100,102) 'rij','g(rij)'
            do i=1,n_bins
                rij=step_bins*(i+0.5_dp)
                ! volumen between bin i+1 and bin i
                volume=(real(i+1,dp)**3-real(i,dp)**3)*(step_bins**3)
                ! number of ideal gas particles in volume
                nid=(1._dp/3._dp)*16._dp*atan(1._dp)*volume*density
                ! normalize g(rij)
                g(i)=g(i)/(real(ngr*n_p,dp)*nid)
                write(100,101) rij,g(i)
            end do
        close(100)
    end subroutine radial_ditribution_function

end module module_mc_lennard_jones