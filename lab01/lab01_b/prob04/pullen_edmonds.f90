program pullen_edmonds
    implicit none
    
    use module_presition
    use module_EDO_segundo_orden
    use module_pullen_edmonds

    n_opt = 4194641 ! aprox 2^22 (n_opt = (8/eps_machine)^(2/5))

    ! RK4_four_eq(RK4_type,n,a1,b1,y1_0,y2_0,y3_0,y4_0,y1_RK4,y2_RK4,y3_RK4,y4_RK4,function_type,input_type)
    call RK4_four_eq(1_sp,n,a1,b1,y1_0,y2_0,y3_0,y4_0,y1_RK4,y2_RK4,y3_RK4,y4_RK4,1_sp,1_sp)

end program pullen_edmonds