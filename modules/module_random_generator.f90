module module_random_generator
    use module_precision
    implicit none

    contains
    ! ========================================================================
    ! PURPOSE: generador ran2
    ! ========================================================================
    function ran0(idum)
        integer(sp), intent(inout)  :: idum
        integer(sp), parameter      :: ia=16807,im=2147483647,iq=127773,ir=2836
        integer(sp)                 :: k
        real(dp),    parameter      :: am=1._dp/real(im,dp)
        real(dp)                    :: ran0
        k=idum/iq
        idum=ia*(idum-k*iq)-ir*k
        if (idum <= 0) idum = idum + im
        ran0 = am*real(idum,dp)
    end function ran0
    ! ========================================================================
    ! PURPOSE: generador ran2
    ! ========================================================================
    function ran2(idum)
        integer(sp) :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
        real(dp)    :: ran2,AM,EPS,RNMX
        parameter (IM1=2147483563_sp,IM2=2147483399_sp,AM=1._dp/IM1,IMM1=IM1-1_sp, &
        IA1=40014_sp,IA2=40692_sp,IQ1=53668_sp,IQ2=52774_sp,IR1=12211_sp,IR2=3791_sp, &
        NTAB=32_sp,NDIV=1_sp+IMM1/NTAB,EPS=1.2e-7_dp,RNMX=1._dp-EPS)
        integer(sp) :: idum2,j,k,iv(NTAB),iy
        save iv,iy,idum2
        data idum2/123456789_sp/, iv/NTAB*0_sp/, iy/0_sp/
        if (idum.le.0) then
          idum=max(-idum,1)
          idum2=idum
          do j=NTAB+8_sp,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
          end do
          iy=iv(1)
        endif
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum.lt.0) idum=idum+IM1
        k=idum2/IQ2
        idum2=IA2*(idum2-k*IQ2)-k*IR2
        if (idum2.lt.0) idum2=idum2+IM2
        j=1+iy/NDIV
        iy=iv(j)-idum2
        iv(j)=idum
        if(iy.lt.1)iy=iy+IMM1
        ran2=min(AM*real(iy,dp),RNMX)
        return
    end function ran2
    ! ========================================================================
    ! PURPOSE
    ! ========================================================================

end module module_random_generator