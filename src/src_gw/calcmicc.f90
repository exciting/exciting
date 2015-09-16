!
    subroutine calcmicc
!
!!DESCRIPTION:
!
!This subroutine calculates the matrix elements $M^i_{cc}(\vec{k},\vec{q})$ 
!(c stands for core-state)
!
!!USES:
    use modinput
    use modmain
    use modgw

!!INPUT PARAMETERS:
    implicit none

!!LOCAL VARIABLES:
    integer(4) :: icg1, ia1, is1, ias1, ic1, l1, m1
    integer(4) :: icg2, ia2, is2, ias2, ic2, l2, m2
    integer(4) :: imix, irm, im, bl, bm
     
    real(8) :: tstart, tend
     
    complex(8) :: angint

!!EXTERNAL ROUTINES: 
    real(8), external :: gaunt

!!REVISION HISTORY:
! 
! Created  Apr 2015 by DIN
!
!EOP
!BOC

    call timesec(tstart)

    !---------------------------------------------------
    ! n => Valence or conduction state, n' => core state
    !---------------------------------------------------
    if (allocated(miccmat)) deallocate(miccmat)
    allocate(miccmat(locmatsiz,ncg,ncg))
    miccmat = zzero

    ! Loop over core states
    do icg1 = 1, ncg
      is1  = corind(icg1,1)
      ia1  = corind(icg1,2)
      ic1  = corind(icg1,3)
      l1   = corind(icg1,4)
      m1   = corind(icg1,5)      
      ias1 = idxas(ia1,is1)
      
      do icg2 = 1, ncg
        is2  = corind(icg2,1)
        ia2  = corind(icg2,2)
        ic2  = corind(icg2,3)
        l2   = corind(icg2,4)
        m2   = corind(icg2,5)      
        ias2 = idxas(ia2,is2)
        
        if (ias1==ias2) then
        
          ! Loop over mixed functions:
          im = 0
          do irm = 1, nmix(ias1)  
            bl = bigl(ias1,irm)
            do bm = -bl, bl
              im = im+1
              imix = locmixind(ias1,im)
              
              if ((iabs(l1-l2)<=bl).and.(l1+l2>=bl)) then
                angint = gaunt(l1,bl,l2,m1,bm,m2)
                if (abs(angint)>1.0d-8) then
                  miccmat(imix,icg1,icg2) = &
                  &  angint*bradketc(ias1,irm,ic1,ic2,1,1)
                end if ! angint
               end if
               
            end do ! bm  
          end do ! irm
          
        end if
      end do ! icg2
      
    enddo ! icg1

    !------------------------------------------------------------------
    ! timing
    call timesec(tend)
    
    !write(*,*) ' miccmat:' 
    !do imix = 1, locmatsiz, locmatsiz/10
    !  do icg1 = 1, ncg
    !    do icg2 = 1, 14
    !      write(*,*) imix, icg1, icg2, minc(imix,icg1,icg2)
    !    end do
    !  end do
    !end do

    return
end subroutine
