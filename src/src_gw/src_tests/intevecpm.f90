!BOP
!
! !ROUTINE: intevecpm
!
! !INTERFACE:
      subroutine intevecpm(ib1,ib2)

! !DESCRIPTION:
!
! This subroutine calculates the integral of the product of two LAPW
! eigenvectors by integration in the MT sphere and by expandding it in the
! mixed basis
!
! !USES:
      
      use modmain
      use modgw

! !INPUT PARAMETERS:

      implicit none

      integer(4), intent(in) :: ib1 ! The band index of the function

      integer(4), intent(in) :: ib2 ! The band index of the function

! !LOCAL VARIABLES:

      integer(4) :: im, irm, l, m, imix
      integer(4) :: ii1, ii2, jatom
      integer(4) :: ia, is, ias
      
      real(8) :: abserr,relerr
      
      real(8) :: epii,intep
      real(8) :: epint(natmtot)
      
      complex(8) :: intmix,emii
      complex(8) :: emint(natmtot)
! 
! !REVISION HISTORY:
!
! Created: 9th. July 2004 by RGA
! Last modified: 15th. July 2004 by RGA
! Revisited: May 2011 by DIN
!
!EOP
!BOC
!------------------
!     MT region
!------------------
      intmix = zzero
      intep = 0.0d0
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
!
!         Integral over products of KS orbitals
!
          read(71,13) ii1, ii2, jatom, epint(jatom)
          intep = intep+epint(jatom)
!          
          if(ii1.ne.ib1)then
            write(*,*)'ii1 = ',ii1,'ib1 = ',ib1
            write(*,*)'ii2 = ',ii2,'ib2 = ',ib2
            write(*,*)'jatom = ',jatom,'latom = ',ias
            stop 'at error ii1.ne.ib1'
          endif  
          if(ii2.ne.ib2)then
            write(*,*)'ii1 = ',ii1,'ib1 = ',ib1
            write(*,*)'ii2 = ',ii2,'ib2 = ',ib2
            write(*,*)'jatom = ',jatom,'latom = ',ias
            stop 'at error ii2.ne.ib2'
          endif  
          if(jatom.ne.ias)then
            write(*,*)'ii1 = ',ii1,'ib1 = ',ib1
            write(*,*)'ii2 = ',ii2,'ib2 = ',ib2
            write(*,*)'jatom = ',jatom,'latom = ',ias
            stop 'error jatom.ne.latom'
          endif  
!
!         Integral over mixed wavefunctions
!
          emint(ias) = zzero
          imix = 0
          do irm = 1, nmix(ias)
            l = bigl(irm,ias)
            do m = -l, l
              imix = imix+1
              im = locmixind(ias,imix)
              emint(ias) = emint(ias)+ &
              &            minmmat(im,ib1,ib2)*conjg(minmmat(im,ib1,ib2))
            end do
          end do 
          intmix = intmix+emint(ias)
!
          abserr = abs(epint(jatom)-real(emint(ias)))
          relerr = abserr/epint(ias)
          write(74,11)ib1,ib2,ias,epint(jatom),real(emint(ias)),abserr,relerr
        
        enddo ! ia
      enddo ! is    

!----------------
!     Interstitial    
!----------------
      read(72,14)ii1,ii2,epii
      if(ii1.ne.ib1)then
        write(*,*)'ii1 = ',ii1,'ib1 = ',ib1
        write(*,*)'ii2 = ',ii2,'ib2 = ',ib2
        stop 'int error ii1.ne.ib1'
      endif  
      if(ii2.ne.ib2)then
        write(*,*)'ii1 = ',ii1,'ib1 = ',ib1
        write(*,*)'ii2 = ',ii2,'ib2 = ',ib2
        stop 'int error ii2.ne.ib2'
      endif  
      intep=intep+epii
!
      emii=zzero
      do im=locmatsiz+1,matsiz
        emii=emii+minmmat(im,ib1,ib2)*conjg(minmmat(im,ib1,ib2))
      enddo  
      intmix=intmix+emii
!
      abserr=abs(epii-real(emii))
      relerr=abserr/epii
      
      write(75,10)ib1,ib2,epii,real(emii),abserr,relerr

!-------------------------------
!     The total error estimate
!-------------------------------
      
      abserr=abs(intep-real(intmix))
      relerr=abserr/intep
      write(76,10)ib1,ib2,intep,real(intmix),abserr,relerr

   10 format(2i4,4d15.7)
   11 format(3i4,4d15.7)

   13 format(3i4,d15.7)
   14 format(2i4,d15.7)
      return
  
      end subroutine
!EOC



