module module_2D_ferromagnetic_ising
    use module_precision;use module_mt19937, only: sgrnd,grnd
    implicit none
    contains

    subroutine initial_spins_configuration(T_init_type,n,aux_matrix_pbc)
        integer(sp), intent(in)    :: T_init_type,n
        integer(sp), intent(inout) :: aux_matrix_pbc(n+2,n+2)

        integer(sp) :: seed,seed_val(8),i,j
        real(dp)    :: nrand

        call date_and_time(values=seed_val)
        seed=seed_val(8)*seed_val(7)*seed_val(6)+seed_val(5);call sgrnd(seed)

        ! creamos matriz de spin con PBC
        select case (T_init_type)
            case(1) ! genero configuracion inicial (random ó T(inicial)=Infinity)
                do j=2,n+1;do i=2,n+1
                    nrand=real(grnd(),dp)
                    if (nrand<0.5_dp) then;aux_matrix_pbc(i,j)=-1_sp
                    else; aux_matrix_pbc(i,j)=1_sp;end if
                end do;end do
                call pbc(aux_matrix_pbc,n,5_sp) ! cambio de bordes
            case(2) ! genero configuracion inicial (ordenada ó T(inicial)=0)
                aux_matrix_pbc(:,:) = 1_sp ! todos los spins "up"
            case(3) ! genero configuracion inicial (ordenada ó T(inicial)=0)
                aux_matrix_pbc(:,:) = -1_sp ! todos los spins "down"
        end select

    end subroutine initial_spins_configuration

    ! Energia en unidades de J_int (coeficiente de interacción >0 => ferromagnetic)
    subroutine average_energy(aux_matrix_pbc,n,U)
        integer(sp), intent(in)    :: n
        integer(sp), intent(in)    :: aux_matrix_pbc(1:n+2,1:n+2)
        real(dp),    intent(inout) :: U
        integer(sp) :: i,j

        U=0._dp
        do j=2,n+1
            do i=2,n+1
                U=U-real(aux_matrix_pbc(i,j)*(aux_matrix_pbc(i,j+1)+aux_matrix_pbc(i,j-1)&
                    & +aux_matrix_pbc(i+1,j)+aux_matrix_pbc(i-1,j)),dp)
            end do
        end do
        U=U*0.5_dp ! suponemos J_int >0 (to compensate over-counting)
    end subroutine average_energy

    ! Delta energético por "flipear" un spin, en unidades de J_int
    function U_delta_adim(spin_value,spin1,spin2,spin3,spin4)
        integer(sp), intent(in) :: spin_value,spin1,spin2,spin3,spin4
        integer(dp)             :: U_delta_adim
        U_delta_adim=2_dp*real(spin_value*(spin1+spin2+spin3+spin4),dp)
    end function U_delta_adim

    ! MAGNETIZACIÓN APROXIMADA (específica ó por unidad de partícula)
    function M_adim(aux_matrix_pbc,n)
        integer(sp), intent(in) :: n,aux_matrix_pbc(n+2,n+2)
        integer(sp)             :: i,j
        real(dp)                :: M_adim
        M_adim=0._dp
        do j=2,n+1;do i=2,n+1
            M_adim=M_adim+real(aux_matrix_pbc(i,j),dp)
        end do;end do
        M_adim=(1._dp/(n*n))*abs(M_adim) ! con valor absoluto
        !M_adim=(1._dp/(n*n))*M_adim ! sin valor absoluto
    end function M_adim

    ! MAGNETIZACIÓN EXACTA
    function M_exact_adim(n,T_adim,Tc_adim)
        integer(sp), intent(in) :: n
        real(dp),    intent(in) :: T_adim,Tc_adim
        real(dp)                :: zz,M_exact_adim
        cond1: if (T_adim<Tc_adim) then
            if (T_adim==0._dp) then;M_exact_adim=1._dp;exit cond1;end if
            if (4*(1._dp/T_adim)<abs(log(tiny(1._dp)))) then
                zz=exp(-4*(1._dp/T_adim)) ! zz=z**2
            else;zz=0._dp;end if
            M_exact_adim=((1+zz)**0.25_dp)*((1-6*zz+zz*zz)**0.125_dp)*(1._dp/sqrt(1._dp-zz))
        else;M_exact_adim=0._dp;end if cond1
    end function M_exact_adim

    ! CALCULO DE RELAJACIÓN CON MÉTODO MONTE CARLO METRÓPOLIS
    subroutine MC_step_relaxation(MC_step_type,n,aux_matrix_pbc,T_adim,U_adim)

        integer(sp), intent(in)    :: MC_step_type,n            ! Tipo de MC relaxation y dimension
        integer(sp), intent(inout) :: aux_matrix_pbc(n+2,n+2)   ! Matriz de spins
        real(dp),    intent(inout) :: U_adim                    ! Energía interna
        real(dp),    intent(in)    :: T_adim                    ! Temperatura de equilibrio
        
        real(dp)    :: deltaU_adim
        real(dp)    :: nrand
        integer(sp) :: seed,seed_val(8),i,j,k

        call date_and_time(values=seed_val)
        seed=seed_val(8)*seed_val(7)*seed_val(6)+seed_val(5);call sgrnd(seed)

        ! hago un paso de Monte Carlo (MC_step)
        select case(MC_step_type)
            case(1) ! exploramos el espacio de configuraciones de forma random
                do k=1,n*n
                    ! floor(x,type) returns the greatest integer less than or equal to X
                    nrand=real(grnd(),dp);i=2_sp+floor(n*nrand,sp) ! random integer from 2 to n+1 (fila)
                    nrand=real(grnd(),dp);j=2_sp+floor(n*nrand,sp) ! random integer from 2 to n+1 (columna)

                    deltaU_adim=U_delta_adim(aux_matrix_pbc(i,j),&
                                            aux_matrix_pbc(i,j+1),aux_matrix_pbc(i,j-1),&
                                            aux_matrix_pbc(i+1,j),aux_matrix_pbc(i-1,j))
                    cond1:if (deltaU_adim<=0._dp) then
                        U_adim=U_adim+deltaU_adim
                        aux_matrix_pbc(i,j)=-aux_matrix_pbc(i,j)
                        call pbc_ij(aux_matrix_pbc,n,i,j) ! actualizamos condiciones de borde
                    else
                        if (T_adim==0._dp) exit cond1
                        nrand=real(grnd(),dp)
                        if (deltaU_adim*(1._dp/T_adim)<abs(log(tiny(1._dp)))) then
                            if (exp(-deltaU_adim*(1._dp/T_adim))>=nrand) then
                                U_adim=U_adim+deltaU_adim
                                aux_matrix_pbc(i,j)=-aux_matrix_pbc(i,j)
                                call pbc_ij(aux_matrix_pbc,n,i,j) ! actualizamos condiciones de borde
                            end if
                        else
                            if (nrand==0._dp) then
                                U_adim=U_adim+deltaU_adim
                                aux_matrix_pbc(i,j)=-aux_matrix_pbc(i,j)
                                call pbc_ij(aux_matrix_pbc,n,i,j) ! actualizamos condiciones de borde
                            end if
                        end if
                    end if cond1
                end do
            case(2) ! exploramos el espacio de configuraciones de forma ordenada
                do j=2,n+1;do i=2,n+1
                    deltaU_adim=U_delta_adim(aux_matrix_pbc(i,j),&
                                            aux_matrix_pbc(i,j+1),aux_matrix_pbc(i,j-1),&
                                            aux_matrix_pbc(i+1,j),aux_matrix_pbc(i-1,j))
                    cond2: if (deltaU_adim<=0._dp) then
                        U_adim=U_adim+deltaU_adim
                        aux_matrix_pbc(i,j)=-aux_matrix_pbc(i,j)
                        call pbc_ij(aux_matrix_pbc,n,i,j) ! actualizamos condiciones de borde
                    else
                        if (T_adim==0._dp) exit cond2
                        nrand=real(grnd(),dp)
                        if (deltaU_adim*(1._dp/T_adim)<abs(log(tiny(1._dp)))) then
                            if (exp(-deltaU_adim*(1._dp/T_adim))>=nrand) then
                                U_adim=U_adim+deltaU_adim
                                aux_matrix_pbc(i,j)=-aux_matrix_pbc(i,j)
                                call pbc_ij(aux_matrix_pbc,n,i,j) ! actualizamos condiciones de borde
                            end if
                        else
                            if (nrand==0._dp) then
                                U_adim=U_adim+deltaU_adim
                                aux_matrix_pbc(i,j)=-aux_matrix_pbc(i,j)
                                call pbc_ij(aux_matrix_pbc,n,i,j) ! actualizamos condiciones de borde
                            end if
                        end if
                    end if cond2
                end do;end do
        end select
    end subroutine MC_step_relaxation

    subroutine pbc(aux_matrix_pbc,n,change_type)
        implicit none
        integer(sp), intent(inout) :: aux_matrix_pbc(n+2,n+2)
        integer(sp), intent(in)    :: n,change_type
        select case(change_type) ! tipo de cambio
            case(1) ! cambiamos primer fila     (borde superior)
                aux_matrix_pbc(1,2:n+1)=aux_matrix_pbc(n+1,2:n+1)
            case(2) ! cambiamos última fila     (borde inferior)
                aux_matrix_pbc(n+2,2:n+1)=aux_matrix_pbc(2,2:n+1)
            case(3) ! cambiamos primer columna  (borde izquierdo)
                aux_matrix_pbc(2:n+1,1)=aux_matrix_pbc(2:n+1,n+1)
            case(4) ! cambiamos última columna  (borde derecho)
                aux_matrix_pbc(2:n+1,n+2)=aux_matrix_pbc(2:n+1,2)
            case(5) ! cambiamos todo            (todos los bordes)
                aux_matrix_pbc(n+2,2:n+1)=aux_matrix_pbc(2,2:n+1)
                aux_matrix_pbc(2:n+1,n+2)=aux_matrix_pbc(2:n+1,2)
                aux_matrix_pbc(1,2:n+1)=aux_matrix_pbc(n+1,2:n+1)
                aux_matrix_pbc(2:n+1,1)=aux_matrix_pbc(2:n+1,n+1)
        end select
    end subroutine pbc

    subroutine pbc_ij(aux_matrix_pbc,n,i,j)
        implicit none
        integer(sp), intent(in)    :: n,i,j
        integer(sp), intent(inout) :: aux_matrix_pbc(n+2,n+2)

        if (i==2.and.j==2) then          ! change right & lower edges
            call pbc(aux_matrix_pbc,n,4_sp);call pbc(aux_matrix_pbc,n,2_sp)
        else if (i==2.and.j==n+1) then   ! change left & lower edges
            call pbc(aux_matrix_pbc,n,3_sp);call pbc(aux_matrix_pbc,n,2_sp)
        else if (i==n+1.and.j==2) then   ! change right & upper edges
            call pbc(aux_matrix_pbc,n,4_sp);call pbc(aux_matrix_pbc,n,1_sp)
        else if (i==n+1.and.j==n+1) then ! change left & upper edges
            call pbc(aux_matrix_pbc,n,3_sp);call pbc(aux_matrix_pbc,n,1_sp)
        end if

        if (i==2.and.j/=2.and.j/=n+1)   call pbc(aux_matrix_pbc,n,2_sp) ! cambio borde inferior
        if (j==2.and.i/=2.and.i/=n+1)   call pbc(aux_matrix_pbc,n,4_sp) ! cambio borde derecho
        if (i==n+1.and.j/=2.and.j/=n+1) call pbc(aux_matrix_pbc,n,1_sp) ! cambio borde superior
        if (j==n+1.and.i/=2.and.i/=n+1) call pbc(aux_matrix_pbc,n,3_sp) ! cambio borde izquierdo

    end subroutine pbc_ij

    ! subrutina para calcular la función de autocorrelación
    subroutine func_autocor(obs,time_index,tau_corr,num_of_terms,&
                        autocor_vector,aux_vector1,mask_vector,aux_vector2,&
                        obs_med,var)
        implicit none

        integer(sp), intent(in)    :: time_index    ! indice del observable
        integer(sp), intent(in)    :: tau_corr      ! tiempo de correlación
        integer(sp), intent(in)    :: num_of_terms  ! menor cant de terminos admisibles (debe ser multiplo de tau_corr)
        real(dp),    intent(in)    :: obs           ! observable
        real(dp),    intent(inout) :: autocor_vector(tau_corr)                      ! vector autocorrelación
        real(dp),    intent(inout) :: aux_vector1(tau_corr),aux_vector2(tau_corr)   ! vectores auxiliares
        real(dp),    intent(inout) :: mask_vector(tau_corr)                         ! vector máscara
        real(dp),    intent(inout) :: obs_med,var                                   ! valor medio y varianza
        integer(sp)                :: i,index                                       ! loop indices
        integer(sp)                :: total_obs_num ! total de observables que deben

        ! calculamos valores primeros y segundos momentos
        obs_med=obs_med+obs ! primeros momentos (acumulación)
        var=var+obs*obs     ! segundos momentos (acumulación)

        ! llenamos por primera vez el vector de dim=tau_corr
        if (time_index<=tau_corr) autocor_vector(time_index)=obs
        if (time_index==tau_corr) then
            aux_vector1(:)=autocor_vector(:)
            aux_vector2(:)=autocor_vector(:)
        end if

        ! determinamos total de observables que van a ingresar
        total_obs_num=num_of_terms+tau_corr

        ! computamos las correlaciones a todo tiempo
        if (time_index>tau_corr.and.time_index<=(total_obs_num)) then
            ! definimos el indice dentro del rango [1,tau_corr]
            if (mod(time_index,tau_corr)==0_sp) then;index=tau_corr
            else;index=mod(time_index,tau_corr);end if
            aux_vector2(:)=cshift(aux_vector2(:),shift=1)
            aux_vector2(tau_corr)=obs
            autocor_vector(index)=autocor_vector(index)+dot_product(aux_vector1(:),aux_vector2(:))
        end if
        if (mod(time_index,tau_corr)==0_sp) aux_vector1(:)=aux_vector2(:)

        ! computamos las correlaciones faltantes
        if (time_index==total_obs_num) then
            mask_vector(:)=1._dp
            aux_vector1(:)=aux_vector2(:)
            do index=1,tau_corr-1
                aux_vector2(:)=cshift(aux_vector2(:),shift=1)
                mask_vector(tau_corr-(index-1))=0._dp
                aux_vector2(:)=aux_vector2(:)*mask_vector(:)
                autocor_vector(index)=autocor_vector(index)+dot_product(aux_vector1(:),aux_vector2(:))
                ! promediamos en el ensamble
                autocor_vector(index)=autocor_vector(index)*(1._dp/real(num_of_terms+tau_corr-index,dp))
            end do
            ! promediamos en el ensamble (tau_corr términos)
            autocor_vector(tau_corr)=autocor_vector(tau_corr)*(1._dp/real(num_of_terms,dp))
            obs_med=obs_med*(1._dp/real(total_obs_num,dp)) ! primeros momentos
            var=var*(1._dp/real(total_obs_num,dp)) ! segundos momentos
            var=var-obs_med*obs_med ! varianza
        end if
    end subroutine func_autocor

end module module_2D_ferromagnetic_ising