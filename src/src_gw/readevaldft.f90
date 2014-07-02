!BOP
! !ROUTINE: readeval
! !INTERFACE:
subroutine readevaldft
! !USES:
    use modmain
    use modgw
    use modmpi

! !DESCRIPTION:
!   Inputs the second-variational eigenvalues and occupancion numbers for
!   valence and core electrons
!
! !REVISION HISTORY:
!   Created May 2006 (RGA)
!   Revisited: 28.04.2011 (DIN)
!EOP
!BOC
    implicit none
! local variables
    integer :: ik,ik0,ist,jst,is,ia,ias,io,iu
    integer :: i,n
    real(8) :: egap, e0, e1
    real(8), parameter :: eps=1.0d-4
    
    integer(4) :: recl
    integer :: nstfv_, nspnfv_
    real(8) :: vkl_(3)
    logical :: exist
    
    character(256) :: filename
    
    if (allocated(evaldft)) deallocate(evaldft)
    allocate(evaldft(nstfv,nkpt))
    evaldft(:,:)=0.d0

!-------------------------------------------------------------------
!   read the KS eigenenergies (only for the irreducible k-points)
!-------------------------------------------------------------------

    filename = "EVALFV_GW.OUT"

!   find the record length
    Inquire (IoLength=Recl) vkl_, nstfv_, nspnfv_
    Do i = 1, 100
      Inquire (File=trim(filename), Exist=Exist)
      If (exist) Then
        Open (70, File=trim(filename), Action='READ', &
        &  Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
        Exit
      Else
        Call system ('sync')
        Write (*,*) '(readevaldft): Waiting for other process to read'
        Call sleep (5)
      End If
    End Do
    Read(70, Rec=1) vkl_, nstfv_, nspnfv_
    close(70)
    
    if (nstfv_.lt.nstfv) then
      write(*,*) 'ERROR(readevaldft) Different number of KS states in EVALFV.OUT!'
      write(*,*) '  nstfv=', nstfv,'  nstfv_=', nstfv_
      stop
    end if

    Inquire (IoLength=Recl) vkl_, nstfv_, nspnfv_, evaldft(:,1)
    Open (70, File=trim(filename), Action='READ', &
   &  Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
    do ik = 1, nkpt
      ik0=idikp(ik)
      Read(70, Rec=ik0) vkl_, nstfv_, nspnfv_, evaldft(:,ik)
    end do ! ik
    Close (70)

!----------------------------------------
! find Fermi energy
!----------------------------------------

    n = int(chgval/2.d0)+10
    call fermi(nkpt,n,evaldft(1:n,:),ntet,tnodes,wtet,tvol, &
    &          chgval,.false.,efermi,egap)
    
!-------------------------------
!   KS band structure analyse
!-------------------------------
    call bandanalysis('KS',ibgw,nbgw,evaldft(ibgw:nbgw,:),efermi)
!
!   Check the range of GW bands [ibgw,nbgw] and 
!   the corresponding number of valence electrons involved in GW calculations
!
    if (ibgw.ge.numin) then
       write(fgw,*) "ERROR(readevaldft): Wrong band interval speciefied"
       write(fgw,*) " ibgw=",ibgw," >= numin=", numin
       stop
    end if   
    if (nbgw.le.nomax) then 
       write(fgw,*)"ERROR(readevaldft): Wrong band interval speciefied"
       write(fgw,*) " nbgw=", nbgw, "<= nomax=", nomax
       stop
    endif 
    nbandsgw=nbgw-ibgw+1
    nvelgw=chgval-2.d0*dble(ibgw-1)
    
    write(fgw,*)
    write(fgw,*)'Maximum number of LAPW states (determined by rgkmax):', nmatmax
    write(fgw,*)'Number of states used in GW:                 '
    write(fgw,*)'    - total                                  ', nstfv
    write(fgw,*)'    - occupied                               ', int(chgval/2.d0)
    write(fgw,*)'    - unoccupied                             ', input%gw%nempty
    
    e0=maxval(evaldft(nstfv,:))
    write(fgw,'("Energy cutoff corresponding to highest unoccupied state [Ha]: ",f12.6)') e0
    write(fgw,'("Energy cutoff corresponding to highest unoccupied state [eV]: ",f12.6)') e0*heV

    write(fgw,*)'Number of gw bands (gw output):              ', nbandsgw
    write(fgw,*)'Range of GW bands:                           ', ibgw, nbgw
    e0=minval(evaldft(ibgw,:))
    e1=maxval(evaldft(nbgw,:))
    write(fgw,'("Range of GW bands [Ha]:    ",2f12.6)') e0, e1
    write(fgw,'("Range of GW bands [eV]:    ",2f12.6)') e0*heV, e1*heV
    write(fgw,*)'Number of valence electrons:                 ', int(chgval) 
    write(fgw,*)'Number of valence electrons treated in GW    ', int(nvelgw) 
!
!   Treatment of the symmetry requires averaging over degenerated states.
!   Array n12dgn contains lower and upper indexes of the degenerated states.
!   This array is used then in the routines to calculate the self-energies.
!
    if(allocated(n12dgn))deallocate(n12dgn)
    allocate(n12dgn(2,nstfv,nkpt))
    
    do ik = 1, nkpt

       if (input%gw%reduceq) then

         ist=1
         do while (ist<=nstfv)
            e0=evaldft(ist,ik)
            ! calculate the number of degenerated states
            n=1; jst=ist+1
            do while (jst<=nstfv)
               if(abs(evaldft(jst,ik)-e0)<eps)then
                  n=n+1
                  jst=jst+1
               else
                  exit
               end if
            end do
            ! indexation of the degenerated states
            do i=0,n-1
               n12dgn(1,ist+i,ik)=ist
               n12dgn(2,ist+i,ik)=ist+n-1
            end do
            ist=ist+n
         enddo ! ist
       
       else
       
         ! no symmetry: degeneracy index is set to 1
         do ist = 1, nstfv
           n12dgn(1,ist,ik)=ist
           n12dgn(2,ist,ik)=ist
         end do
       
       end if ! reduceq
       
       ! Debug information
       if (debug) then

         write(55,*)
         write(55,*)'Degeneracy KS bands: ik = ', ik
         do ist=1,nstfv
            write(55,*)'ist, n12dgn: ', ist, n12dgn(:,ist,ik)
         end do
         write(55,*)
         
       end if ! debug
    
    enddo ! ik

return
end subroutine
!EOC
