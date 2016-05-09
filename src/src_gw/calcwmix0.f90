!BOP
!
!!ROUTINE: calcwmix0
!
!!INTERFACE:
!
subroutine calcwmix0
!
!!DESCRIPTION:
!
!This subroutine calculate the matrix elements $\mathcal{W}^i_0$
!      
!!USES:
    use modinput
    use modmain,               only : pi, zzero, idxas
    use modgw,                 only : Gqset, fdebug
    use mod_product_basis,     only : locmatsiz, matsiz, mpwipw, mbindex, rtl
    use mod_coulomb_potential, only : wi0
    use mod_misc_gw,           only : vi

!!LOCAL VARIABLES:
    implicit none
    integer(4) :: imix, l1, irm
    integer(4) :: ia, is, ias
    integer(4) :: igq
    real(8)    :: fact
      
!!REVISION HISTORY:
!
! Created 11.02.05 by RGA
! Revisited June 2011 by DIN
!
!EOP
!BOC
    if (allocated(wi0)) deallocate(wi0)
    allocate(wi0(matsiz))
    wi0(1:matsiz) = zzero

    !---------
    ! MT part
    !---------
    ! Y_{00}=1/sqrt(4*pi)
    fact = dsqrt(4.0d0*pi*vi)
    do imix = 1, locmatsiz
      l1  = mbindex(imix,4)
      if (l1==0) then
        is  = mbindex(imix,1)
        ia  = mbindex(imix,2)
        irm = mbindex(imix,3)
        ias = idxas(ia,is)
        wi0(imix) = fact*rtl(irm,ias)
      end if
    end do ! imix

    !---------
    ! PW part 
    !---------
    do igq = 1, Gqset%ngk(1,1)
      imix = locmatsiz+igq
      wi0(imix) = mpwipw(igq,1)
    end do

    if (input%gw%debug) then
      write(fdebug,*) '---- calcwmix0 < mixind    wi0 >'
      do imix = 1, matsiz, 10
        write(fdebug,*) imix, wi0(imix)
      end do
    end if

end subroutine
!EOC                
