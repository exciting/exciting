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
      use modmpi
            
! !LOCAL VARIABLES:

      implicit none

      integer(4) :: iq, iqp
      integer(4) :: ik, ikp,ikpqp
      integer(4) :: recl
      character(128)::sbuffer
      real(8)    :: tq1, tq2, tq11, tq22
      complex(8),allocatable::buffer(:)
	  integer :: COMM_LEVEL2
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
   
    call barrier()
 
!=========================================================================================
!     Calculate the dielectric function matrix
!=========================================================================================
    
      if (.not.input%gw%reps) then
      
        call cpu_time(tq1)
       
        if(allocated(epsilon))deallocate(epsilon)
        allocate(epsilon(matsizmax,matsizmax,nomeg))
	    epsilon=zzero
        if(allocated(inveps))deallocate(inveps)
        allocate(inveps(matsizmax,matsizmax,nomeg))
		inveps=zzero
!--------------------------------------------------------------------------
!       The direct access file to store the values of the inverse dielectric matrix
!--------------------------------------------------------------------------
        recl=16*(matsizmax*matsizmax*nomeg)
       !this buffer is used only for making INVEPS files more comparable
        allocate(buffer(matsizmax*matsizmax*nomeg))
        buffer=zzero
        !each process opens its own file except there is no q-point left
        write(sbuffer,*)rank
       if (rank.lt.nqpt) open(44,file='INVEPS'//trim(adjustl(sbuffer))//'.OUT',action='WRITE',form='UNFORMATTED', &
       &  access='DIRECT',status='REPLACE',recl=recl)
!
!       Loop over q-points
!        
        call boxmsg(fgw,'-','Calculation of the dielectric matrix')
        call barrier

         write(*,*)"On proc",rank,"do qp",firstofset(mod(rank,nqpt),nqpt)," to ", lastofset(mod(rank,nqpt),nqpt)
#ifdef MPI
		! create  communicator object for each qpoint
         call MPI_COMM_SPLIT(MPI_COMM_WORLD,firstofset(mod(rank,nqpt),nqpt),  rank/nqpt, COMM_LEVEL2,ierr)
        ! write(*,*) "proc ",rank,"does color ", firstofset(mod(rank,nqpt),nqpt), "with key(rank)", rank/nqpt
#endif
		!each process does a subset
        do iqp = firstofset(mod(rank,nqpt),nqpt), lastofset(mod(rank,nqpt),nqpt)


          call cpu_time(tq11)
!      
!         Calculate the interstitial mixed basis functions
!  
          iq=idikpq(iqp,1)
          Gamma=gammapoint(iq)
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
          call setbarcev(barcevtol)



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
          call calcepsilon(iqp,COMM_LEVEL2)
!
!         Calculate inverse of the dielectric function
!   
          call calcinveps(iqp)
          
          ! store the inverse dielmat into file INVEPS.OUT
          if (rank.lt.nqpt) then
          	 write(44,rec=iqp-firstofset(rank,nqpt)+1) buffer !fill file with zeroes
        	 write(44,rec=iqp-firstofset(rank,nqpt)+1) inveps
        	! write(1100+iqp,*)epsilon
        	! write(900+iqp,*)inveps
   			 write(fgw,*)"write iqp ",iqp,"in file: ", rank,"at rec:" ,iqp-firstofset(rank,nqpt)+1
		  endif
          if((iqp.eq.1) )then
			if( rank.eq.0)then
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
			endif
            deallocate(head)
            deallocate(epsw1,epsw2)
            deallocate(emac)
          endif

          call cpu_time(tq22)
          
          write(fgw,*) 'q-point =', iqp, '    takes ', tq22-tq11, ' CPU seconds'

        end do ! iqp

        if (rank.lt.nqpt) close(44)
        call barrier

        deallocate(epsilon)
        deallocate(inveps)

        
        
        call cpu_time(tq2)
        write(fgw,*)
        call write_cputime(fgw,tq2-tq1,'TOTAL CPU TIME:')
        write(fgw,*)
        
      end if
      
      !debug info
      if(1.eq.0)then
      do iqp=1,nkpt
      call getinveps(iqp)
      write(8888,*)inveps
      end do
      end if

!=========================================================================================
!                                  MAIN LOOP
!=========================================================================================

!     The file to store intermediate (q-dependent) values of \Sigma_{x,c}
      open(96,file='ADDSELFE'//trim(adjustl(sbuffer))//'.OUT',form='FORMATTED',status='UNKNOWN')

!     Calculate the integrals to treat the singularities at G+q->0
      call setsingc
      call barrier
  	  ikpqp=0 ! initialize counter for k-points
      do ikp = 1, nkpt ! IBZ

         do iqp = 1, nkptq(ikp)  ! EIBZ(k)
         ikpqp=ikpqp+1
		! decide if point is done by this process
         if (mod(ikpqp-1,procs).eq.rank) then
         write(*,*) "do q",iqp,"and k", ikp,"as",ikpqp, "on",rank
           
           call cpu_time(tq1)
           
           ik=idikp(ikp)
           iq=idikpq(iqp,ikp)
           Gamma=gammapoint(iq)
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
           call setbarcev(barcevtol)
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


        endif
        end do ! iqp
      
      end do ! ikp
!##############################      
#ifdef MPI
 call MPI_ALLREDUCE(MPI_IN_PLACE, selfec, (nbgw-ibgw+1) *nkpt*nomeg,  MPI_DOUBLE_COMPLEX,  MPI_SUM,&
       & MPI_COMM_WORLD, ierr)
 call MPI_ALLREDUCE(MPI_IN_PLACE, selfex, (nbgw-ibgw+1) *nkpt,  MPI_DOUBLE_COMPLEX,  MPI_SUM,&
       & MPI_COMM_WORLD, ierr)
#endif
      if (allocated(minmmat))deallocate(minmmat)
      if(iopcore.eq.0 .and. allocated(mincmat))deallocate(mincmat)
      close(96)
!
!     Write the exchange term to file
!      
	if (rank.eq.0) then
      open(92,file='SELFX.OUT',form='UNFORMATTED',status='UNKNOWN')
      write(92) ibgw, nbgw, nkpt, selfex
      close(92)
!
!     Write the correlation term to file
!      
      open(93,file='SELFC.OUT',form='UNFORMATTED',status='UNKNOWN')
      write(93) ibgw, nbgw, nkpt, nomeg, selfec
      close(93)
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
      

      close(52) ! STRCONST.OUT
    endif  
        close(46) ! ADDSELFE.OUT
      return
      end subroutine gwcycle
!EOC      
