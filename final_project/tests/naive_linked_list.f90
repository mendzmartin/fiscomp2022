! gfortran -o naive_linked_list.o ../../modules/module_precision.f90 naive_linked_list.f90 && ./naive_linked_list.o
program naive_linked_list
    use module_precision
    implicit none

    use module_precision;use module_md_lennard_jones
    implicit none
    integer(sp), parameter   :: n_p=256_sp                             ! cantidad de partículasa
    real(dp),    parameter   :: delta_time=0.005_dp                    ! paso temporal
    integer(sp), parameter   :: time_eq=1000_sp,time_scal=50_sp,&      ! pasos de equilibración y de escaleo de veloc.
                                time_run=1000_sp                       ! pasos de evolucion en el estado estacionario
    real(dp),    parameter   :: T_adim_ref=1.1_dp                      ! temperatura de referencia adimensional
    real(dp),    parameter   :: density=0.8_dp                         ! densidad (particulas/volumen)
    real(dp),    parameter   :: r_cutoff=2.5_dp,mass=1._dp             ! radio de corte de interacciones y masa     
    real(dp),    allocatable :: x_vector(:),y_vector(:),z_vector(:)    ! componentes de las posiciones/particula
    real(dp),    allocatable :: vx_vector(:),vy_vector(:),vz_vector(:) ! componentes de la velocidad/particula
    real(dp),    allocatable :: force_x(:),force_y(:),force_z(:)       ! componentes de la fuerza/particula
    integer(sp) :: imap,index_x,index_y,index_z,icell
    integer(sp) :: i,j
    integer(sp),allocatable :: head(m*m*m),list(n_p)

    m=
    allocate(head(m*m*m),list(n_p))

    do icell=1,m*m*m ! celdas
        i = head(icell)
        do while (i/=0) !
            j = list(i)
            do while (j/=0) ! en la celda
                . . . calculo para el par i,j
                . . . por ej fxi = fxi + . . .
                . . . por ej fx(j)= fx(j)+ . . .
                j = list(j)
            end do
            jcell0 = 13*(icell - 1)
            do nabor=1,13 ! celdas vecinas
                jcell = MAP(jcell0 + nabor)
                j = head(jcell)
                do while (j/=0) ! en la celda
                    . . . calculo para el par i,j
                    . . . por ej fxi = fxi + . . .
                    . . . por ej fx(j)= fx(j)+ . . .
                    j = list(j)
                end do
            end do
            . . . por ej fx(i) = fxi
            i = list(i)
        end do
    end do ! celdas       

end program

subroutine maps(m,map)
    use module_precision
    implicit none
    integer(sp), intent(in) :: m
    integer(sp), intent(inout) :: map(13*m*m*m)
    integer(sp) :: index_x,index_y,index_z,imap,icell

    do index_x=1,m
        do index_y=1,m
            do index_z=1,m
                imap=(icell(index_x,index_y,index_z,m)-1)*13
                map(imap+1)=icell(index_x+1,index_y,index_z,M)
                map(imap+2)=icell(index_x+1,index_y+1,index_z,M)
                map(imap+3)=icell(index_x,index_y+1,index_z,M)
                map(imap+4)=icell(index_x-1,index_y+1, index_z,M)
                map(imap+5)=icell(index_x+1,index_y,index_z-1,M)
                map(imap+6)=icell(index_x+1,index_y+1,index_z-1,M)
                map(imap+7)=icell(index_x,index_y+1,index_z-1,M)
                map(imap+8)=icell(index_x-1,index_y+1,index_z-1,M)
                map(imap+9)=icell(index_x+1,index_y,index_z+1,M)
                map(imap+10)=icell(index_x+1,index_y+1,index_z+1,M)
                map(imap+11)=icell(index_x,index_y+1,index_z+1,M)
                map(imap+12)=icell(index_x-1,index_y+1,index_z+1,M)
                map(imap+13)=icell(index_x,index_y,index_z+1,M)
            end do
        end do
    end do
end subroutine maps

function icell(index_x,index_y,index_z,M)
    use module_precision
    implicit none
    integer(sp), intent(in) :: index_x,index_y,index_z,M
    integer(sp) :: icell
    icell=1+mod(index_x-1+M,M)+mod(index_y-1+M,M)*M+mod(index_z-1+M,M)*M*M
end function icell

subroutine links(n_p,m,L,head,list,x_vector,y_vector,z_vector)
    use module_precision
    implicit none
    integer(sp), intent(in)    :: m,n_p
    real(dp),    intent(in)    :: L,x_vector(n_p),y_vector(n_p),z_vector(n_p)
    integer(sp), intent(inout) :: head(m*m*m),list(n_p)
    real(dp)                   :: Lc_inv  ! inversa del lado de la celda
    integer(sp)                :: i,icell 
    ! celda de -L/2 a L/2
    head(:)=0_sp
    Lc_inv=real(m,dp)*(1._dp/L)
    do i=1,n_p
        icell=1+int((x_vector(i)+0.5_dp*L)*Lc_inv,sp)+&
            int((y_vector(i)+0.5_dp*L)*Lc_inv,sp)*m+&
            int((z_vector(i)+0.5_dp*L),sp)*m*m
        list(i)=head(icell)
        head(icell)=i
    end do
end subroutine links

