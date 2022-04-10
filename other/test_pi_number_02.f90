program print_pi
! use iso_fortran_env, sp=>real32, dp=>real64, qp=>real128

  integer, parameter :: sp = selected_real_kind(6, 37)
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer, parameter :: qp = selected_real_kind(33, 4931)

  real(kind=sp), parameter :: pi_sp = 3.14159265358979323846264338327950288_sp
  real(kind=dp), parameter :: pi_dp = 3.14159265358979323846264338327950288_dp
  real(kind=qp), parameter :: pi_qp = 3.14159265358979323846264338327950288_qp
  
  write(*,'("SP "A17)') "3.14159265358..."
  write(*,'(F17.11)') pi_sp
  write(*,'(F17.11)')        acos(-1.0_sp)
  write(*,'(F17.11)') 2.0_sp*asin( 1.0_sp)
  write(*,'(F17.11)') 4.0_sp*atan2(1.0_sp,1.0_sp)
  write(*,'(F17.11)') 3.0_sp*acos(0.5_sp)
  write(*,'(F17.11)') 6.0_sp*asin(0.5_sp)

  write(*,'("DP "A26)') "3.14159265358979323846..."
  write(*,'(F26.20)') pi_dp
  write(*,'(F26.20)')        acos(-1.0_dp)
  write(*,'(F26.20)') 2.0_dp*asin( 1.0_dp)
  write(*,'(F26.20)') 4.0_dp*atan2(1.0_dp,1.0_dp)
  write(*,'(F26.20)') 3.0_dp*acos(0.5_dp)
  write(*,'(F26.20)') 6.0_dp*asin(0.5_dp)

  write(*,'("QP "A44)') "3.14159265358979323846264338327950288419..."
  write(*,'(F44.38)') pi_qp
  write(*,'(F44.38)')        acos(-1.0_qp)
  write(*,'(F44.38)') 2.0_qp*asin( 1.0_qp)
  write(*,'(F44.38)') 4.0_qp*atan2(1.0_qp,1.0_qp)
  write(*,'(F44.38)') 3.0_qp*acos(0.5_qp)
  write(*,'(F44.38)') 6.0_qp*asin(0.5_qp)

  write(*,'(F17.11)') sin(pi_sp)
  write(*,'(F26.20)') sin(pi_dp)
  write(*,'(F44.38)') sin(pi_qp)

end program print_pi
