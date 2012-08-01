!BOP
!
! !ROUTINE: excycle
!
! !INTERFACE:
      subroutine excycle
      
! !DESCRIPTION:
!
! This subroutine calculates the exchange self-energy
!
! !USES:

      use modmain
      use modgw
            
! !LOCAL VARIABLES:

      implicit none

      integer(4) :: iq, iqp, ik, ikp
      real(8)    :: tq1,tq2
      
! !REVISION HISTORY:
!
! Created 16.09.2005 by RGA
! Revisited: June 2011 by DIN
! Symmetry added: Dec 2011 by DIN
!      
!EOP
!BOC     
      call boxmsg(fgw,'#','START EXCYCLE')
!
!     Allocate arrays needed for the calculation of the self-energy
!
      allocate(selfex(ibgw:nbgw,nkpt))
      selfex(:,:)=zzero
      
!     Calculate the integrals to treat the singularities at G+q->0
      call setsingc
!
!     The file to store intermediate (q-dependent) values of \Sigma_{x,c}
!
      open(96,file='ADDSELFE.OUT',form='FORMATTED',status='UNKNOWN')

!------------------------------------------------------
!     Main loop
!------------------------------------------------------  
      do ikp = 1, nkpt ! IBZ
      
         ik=idikp(ikp)
         
         do iqp = 1, nkptq(ikp)  ! EIBZ(k)

           call cpu_time(tq1)
           
           iq=idikpq(iqp,ikp)
!
!          Set the size of the basis for the corresponding q-point
!        
           matsiz=locmatsiz+ngq(iq)
           write(fgw,101) ik, iq, locmatsiz, ngq(iq), matsiz
!      
!          Calculate the interstitial mixed basis functions
!       
           call diagsgi(iq)
!      
!          Calculate the transformation matrix between pw's and the mixed basis functions
!
           call calcmpwipw(iq)
! 
!          Calculate the bare coulomb matrix
!
           call calcbarcmb(iq)

           if( barcevtol.gt.0.d0) then
             write(6,'(a,f8.4)')" -Use reduced basis: barcevtol=",&
          &                        barcevtol 
           endif
           call setbarcev(iq,barcevtol)
!        
!          Calculate the Minm(k,q) matrix elements for given k and q
!        
           call expand_prods(ik,iq,1)
!
!          Calculate the exchange self-energy
!
           call calcselfx(ikp,iqp)
        
           call cpu_time(tq2)
           write(fgw,*)
           call write_cputime(fgw,tq2-tq1,'Q-POINT CYCLE TAKES')
           write(fgw,*)
  
        end do ! iqp

      end do ! ikp
      
      deallocate(minmmat)
      if(wcore)deallocate(mincmat)
      close(96)
!
!     Write the exchange term to file
!      
      call writeselfx
!
!     Calculate the diagonal matrix elements of the DFT exchange-correlation potential
!
      call calcvxcnn
!
!     Calculate the quasiparticle energies
!      
      call calceqpx
      
  101 format(10x,/, &
     &       'Data for (k,q)-point:',2i4,//,10x,'Mixed basis:',/,10x, &
     &       'Number of atomic basis functions:       ',i4,/,10x,     &
     &       'Number of interstitial basis functions: ',i4,/,10x,     &
     &       'Total number of basis functions:        ',i4,/)
      
      return
      end subroutine excycle
!EOC      
