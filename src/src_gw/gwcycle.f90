!BOP
!
! !ROUTINE: gwcycle
!
! !INTERFACE:
      
      subroutine gwcycle
      
! !DESCRIPTION:
!
! This subroutine performs one gw cycle and calculates the corresponding
! quasiparticle energies.
!
! !USES:

      use modmain
      use modgw
            
! !LOCAL VARIABLES:

      implicit none

      integer(4) :: iq, iqp
      integer(4) :: ik, ikp
      integer(4) :: recl
      
      real(8)    :: tq1, tq2, tq11, tq22

! !REVISION HISTORY:
!
! Created 16.09.2005 by RGA
! Revisited: June 2011 by DIN
! Symmetry added: Dec 2011 by DIN
!
!EOP
!BOC     
      
      call boxmsg(fgw,'#','START GWCYCLE')
!
!     allocate the arrays needed for the calculation of
!     the self energy
!
      allocate(selfex(ibgw:nbgw,nkpt))
      allocate(selfec(ibgw:nbgw,nkpt,nomeg))
      selfex(:,:)=zzero
      selfec(:,:,:)=zzero
!
!     Calculate the integration weights using the linearized tetrahedron method
!
      call kintw

!=========================================================================================
!     Calculate the momentum matrix elements
!=========================================================================================
      if(.not.input%gw%rpmat)then
        call calcpmat        ! <--- original (RGA's) version
        !call calcpmatgw     ! <--- modified xs (SAG's) version
      else
        write(fgw,*)'PMAT and PMATC are read from file'
      end if

!=========================================================================================
!     Calculate the dielectric function matrix
!=========================================================================================
      if (.not.input%gw%reps) then
      
        call cpu_time(tq1)
       
        if(allocated(epsilon))deallocate(epsilon)
        allocate(epsilon(matsizmax,matsizmax,nomeg))

        if(allocated(inveps))deallocate(inveps)
        allocate(inveps(matsizmax,matsizmax,nomeg))

!--------------------------------------------------------------------------
!       The direct access file to store the values of the inverse dielectric matrix
!--------------------------------------------------------------------------
        recl=16*(matsizmax*matsizmax*nomeg)
        open(44,file='INVEPS.OUT',action='WRITE',form='UNFORMATTED', &
       &  access='DIRECT',status='REPLACE',recl=recl)
!
!       Loop over q-points
!        
        call boxmsg(fgw,'-','Calculation of the dielectric matrix')

        do iqp = 1, nqpt
        
          call cpu_time(tq11)
!      
!         Calculate the interstitial mixed basis functions
!
          iq=idikpq(iqp,1)
          matsiz=locmatsiz+ngq(iq)
!
!         Interstitial mixed basis functions
!
          call diagsgi(iq)
          call calcmpwipw(iq)
!
!         Calculate the bare coulomb potential matrix and its square root
!
          call calcbarcmb(iq)
!
!         Reduce the basis size by choosing eigenvectors of barc with 
!         eigenvalues larger than evtol
!
          call setbarcev(iq,barcevtol)
!
!         Calculate the q-dependent integration weights
!   
          if(convflg.eq.0)then
            call qdepwsum(iq)
          else
            call qdepwtet(iq)
          endif  
!
!         Calculate the dielectric finction
!   
          call calcepsilon(iqp)
!
!         Calculate inverse of the dielectric finction
!   
          call calcinveps(iqp)
          
          ! store the inverse dielmat into file INVEPS.OUT
          write(44,rec=iqp) inveps

          if(iqp.eq.1)then
            ! store the data used for treating q->0 singularities
            open(42,file='INVHEAD.OUT',action='WRITE',form='UNFORMATTED',status='REPLACE')
            write(42)head
            close(42)
            open(43,file='INVWING1.OUT',action='WRITE',form='UNFORMATTED',status='REPLACE')
            write(43)epsw1
            close(43)
            open(43,file='INVWING2.OUT',action='WRITE',form='UNFORMATTED',status='REPLACE')
            write(43)epsw2
            close(43)
            deallocate(head)
            deallocate(epsw1,epsw2)
            deallocate(emac)
          endif
          
          call cpu_time(tq22)
          
          write(fgw,*) 'q-point =', iqp, '    takes ', tq22-tq11, ' CPU seconds'

        end do ! iqp
        
        deallocate(epsilon)
        deallocate(inveps)

        close(44)
        
        call cpu_time(tq2)
        write(fgw,*)
        call write_cputime(fgw,tq2-tq1,'TOTAL CPU TIME:')
        write(fgw,*)
        
      end if
      
!=========================================================================================
!                                  MAIN LOOP
!=========================================================================================

!     The file to store intermediate (q-dependent) values of \Sigma_{x,c}
      open(96,file='ADDSELFE.OUT',form='FORMATTED',status='UNKNOWN')

!     Calculate the integrals to treat the singularities at G+q->0
      call setsingc

      do ikp = 1, nkpt ! IBZ
         
         do iqp = 1, nkptq(ikp)  ! EIBZ(k)
         
           call cpu_time(tq1)
           
           ik=idikp(ikp)
           iq=idikpq(iqp,ikp)
!
!          Set the size of the basis for the corresponding q-point
!        
           matsiz=locmatsiz+ngq(iq)
           write(fgw,101) ikp, iq, locmatsiz, ngq(iq), matsiz
!      
!          Calculate the interstitial mixed basis functions
!       
           call diagsgi(iq)
!      
!          Calculate the transformation matrix between pw's and the mixed basis functions
!
           call calcmpwipw(iq)
! 
!          Calculate the bare coulomb potential matrix
!
           call calcbarcmb(iq)
           call setbarcev(iq,barcevtol)
!        
!          Calculate the Minm(k,q) matrix elements for given k and q
!        
           call expand_prods(ik,iq,1)
!
!          Add the corresponding q-term to the exchange term of the self-energy
!
           call calcselfx(ikp,iqp)
!
!          Calculate the dielectric function
!
           call getinveps(iq)
!
!          Add the corresponding q-term to the correlation term of the self-energy
!
           call calcselfc(ikp,iqp)
        
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
!     Write the correlation term to file
!      
      call writeselfc
!
!     Calculate the diagonal matrix elements of the DFT exchange-correlation potential
!
      call calcvxcnn
!
!     Calculate the quasiparticle energies
!      
      call calceqp

  101 format(10x,/, &
     &       'Data for (k,q)-point:',2i4,//,10x,'Mixed basis:',/,10x, &
     &       'Number of atomic basis functions:       ',i4,/,10x,     &
     &       'Number of interstitial basis functions: ',i4,/,10x,     &
     &       'Total number of basis functions:        ',i4,/)
      
      close(46) ! ADDSELFE.OUT
      close(52) ! STRCONST.OUT
      
      return
      end subroutine gwcycle
!EOC      
