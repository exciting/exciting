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
      use modmpi
            
! !LOCAL VARIABLES:

      implicit none

      integer(4) :: iq, iqp, ik, ikp
      real(8)    :: tq1,tq2
      character(128) :: sbuffer
      complex(8),allocatable :: buffer(:)
      integer :: COMM_LEVEL2
      integer :: ikpqp
      
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

!     The file to store intermediate (q-dependent) values of \Sigma_{x,c}
#ifdef MPI
	  write(sbuffer,*)rank
      open(96,file='ADDSELFE'//trim(adjustl(sbuffer))//'.OUT',form='FORMATTED',status='UNKNOWN')
      if (rank.eq.0 .and. input%gw%debug) write(*,*)"Calculate Self Energy"
#endif
#ifndef MPI
      open(96,file='ADDSELFE.OUT',form='FORMATTED',status='UNKNOWN')
#endif

!     Calculate the integration weights using the linearized tetrahedron method
      call kintw

!     Calculate the integrals to treat the singularities at G+q->0
      call setsingc
      
      call barrier

!------------------------------------------------------
!     Main loop
!------------------------------------------------------
      ikpqp=0  
      do ikp = 1, nkpt ! IBZ
      
         
         do iqp = 1, nkptq(ikp)  ! EIBZ(k)
           
           ikpqp=ikpqp+1
	   ! decide if point is done by this process
           if (mod(ikpqp-1,procs).eq.rank) then
#ifdef MPI
             if (input%gw%debug)write(*,*) "do q",iqp,"and k", ikp,"as number",ikpqp, "on",rank
#endif
           
             call cpu_time(tq1)
           
             ik=idikp(ikp)
             iq=idikpq(iqp,ikp)
             Gamma=gammapoint(iq)
!
!            Set the size of the basis for the corresponding q-point
!        
             matsiz=locmatsiz+ngq(iq)
             write(fgw,101) ik, iq, locmatsiz, ngq(iq), matsiz
!      
!            Calculate the interstitial mixed basis functions
!       
             call diagsgi(iq)
!      
!            Calculate the transformation matrix between pw's and the mixed basis functions
!
             call calcmpwipw(iq)
! 
!            Calculate the bare coulomb matrix
!
             call calcbarcmb(iq)
!        
!            Calculate the Minm(k,q) matrix elements for given k and q
!        
             call expand_prods(ik,iq,1)
!
!            Calculate the exchange self-energy
!
             call calcselfx(ikp,iqp)
        
             call cpu_time(tq2)
             write(fgw,*)
             call write_cputime(fgw,tq2-tq1,'Q-POINT CYCLE TAKES')
             write(fgw,*)
           
           end if
  
        end do ! iqp

      end do ! ikp

!##############################      
#ifdef MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE, selfex, (nbgw-ibgw+1)*nkpt, &
      &  MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
      
      deallocate(minmmat)
      if(iopcore.eq.0)deallocate(mincmat)
      close(96)
!
!     Write the exchange term to file
!      
      if (rank.eq.0) then
      open(92,file='SELFX.OUT',form='UNFORMATTED',status='UNKNOWN')
      write(92) ibgw, nbgw, nkpt, selfex
      close(92)
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
      close(52) ! STRCONST.OUT
      endif  
      close(46) ! ADDSELFE.OUT
#ifdef MPI
      call  cat_logfiles('ADDSELFE')
#endif
      
      return
      end subroutine excycle
!EOC      
