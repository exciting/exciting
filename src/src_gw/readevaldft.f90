!BOP
! !ROUTINE: readeval
! !INTERFACE:
subroutine readevaldft
! !USES:
    use modmain
    use modgw

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
    integer :: ik,ist,jst,is,ia,ias,io,iu
    integer :: i,n
    integer :: ik_, ist_, nkpt_, nstsv_
    real(8) :: egap
    real(8) :: t1, vkl_(3), e0
    real(8), parameter :: eps=1.0d-4
    
    real(8), external  :: idos, dostet
    external fermi

!   allocate second-variational arrays
    if (allocated(evaldft)) deallocate(evaldft)
    allocate(evaldft(nstsv,nkpt))

!   read the valence eigenvalues
    open(50, File='EIGVAL'//trim(filext), Action='READ', Form='FORMATTED')
    read(50,'(I6)') nkpt_
    if (nkpt_.ne.nkpt) stop "readeval: EIGVAL.OUT has different number of k-points"
    read(50,'(I6)') nstsv_
    if (nstsv_.ne.nstsv) stop "readeval: EIGVAL.OUT has different number of states"
    do ik = 1, nkpt
       read(50,*)
       read(50,'(I6,3G18.10)') ik_, vkl_(:)
       t1 = Abs(vkl(1,ik)-vkl_(1)) + Abs(vkl(2, ik)-vkl_(2)) + &
   &        Abs(vkl(3, ik)-vkl_(3))
       if (t1.gt.eps) then
         write (*,*)
         write (*, '("Error(readevaldft): differing vectors for k-point ",I8)') ik
         write (*, '(" current    : ",3G18.10)') vkl (:, ik)
         write (*, '(" EIGVAL.OUT : ",3G18.10)') vkl_
         write (*,*)
         stop
       end if
       read(50,*)
       do ist = 1, nstsv
          read(50,'(I6,2G18.10)') ist_, evaldft(ist,ik), occsv(ist,ik)
       end do
       read(50,*)
    end do
    close (50)

!   read Fermi energy from file    
!    call readfermi
!    write(fgw,*)'(readevaldft) Fermi energy EXCITING [Ha]: ', efermi

!   The fermi energy calculated using LIBBZINT
    call fermi(nkpt,nstfv,evaldft,nirtet,tndi,wirtet,tvol, &
               chgval,.false.,efermi,egap)
    write(fgw,*)
    write(fgw,*)'(readevaldft)     Fermi energy (KS) [Ha]: ', efermi
    write(fgw,*)'(readevaldft)         Band gap (KS) [Ha]: ', egap
    write(fgw,*)

!   KS band structure analyse
    call bandanaly(1,nstfv,nkpt,vkl,evaldft,efermi,"KS",fgw)

!----------------------------------------
! shift all energies to make efermi=0 
!----------------------------------------
    do ik=1,nkpt
      do ist=1,nstfv
        evaldft(ist,ik)=evaldft(ist,ik)-efermi
      enddo
    enddo     

! Core energies are already generated by gencore.f90
! Another option would be to read the data from EVALCORE.OUT
    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        do ist=1,ncore(is)
          evalcr(ist,ias)=evalcr(ist,ias)-efermi
        end do
      end do
    end do

    efermi=0.0d0

!----------------------------------------
! determine additional parameters
!----------------------------------------
    
!   HOMO, LUMO band indices
    maxoccband=0
    minunoband=1000 
    do ik = 1, nkpt
       io=0
       iu=1000
       do ist=1,nstfv
          if(evaldft(ist,ik).lt.efermi)then
             if(ist.gt.io)io=ist
          else
             if(ist.lt.iu)iu=ist
          endif    
       enddo ! ist
       if(io.gt.maxoccband)maxoccband=io
       if(iu.lt.minunoband)minunoband=iu
    enddo ! ik   

    if(maxoccband.lt.1)then
       write(*,'("Error(readeval): maxoccband < 1 : ", i8)')maxoccband
       write(*,*)
       stop
    elseif(maxoccband.ge.nstfv)then
       write(*,'("Error(readeval): maxoccband > nstfv : ", i8)')maxoccband
       write(*,*)
       stop
    endif
    if(minunoband.lt.1)then
       write(*,'("Error(readeval): minunoband < 1 : ", i8)')minunoband
       write(*,*)
       stop
    elseif(minunoband.ge.nstfv)then
       write(*,'("Error(readeval): minunoband > nstfv : ", i8)')minunoband
       write(*,*)
       stop
    endif
    write(fgw,'("Highest occupied band:  ",i5)')maxoccband
    write(fgw,'("Lowest unoccupied band: ",i5)')minunoband 

!
!   Determine the range of GW bands ibgw..nbgw and 
!   the number of valence electrons included in gw bands
!
    if(ibgw.le.0) ibgw=1
    if(nbgw.le.0) nbgw=nstfv
    if(ibgw.ge.minunoband .or. nbgw.le.maxoccband) then 
       write(6,*) "(readevaldft) WARNING: wrong range of gw bands!!!"
       write(6,*) " ibgw >= minunoband=",ibgw,minunoband
       write(6,*) " or"
       write(6,*) " nbgw <= maxoccband=",nbgw,maxoccband
       ibgw=1
       nbgw=nstfv
    endif 
    nbandsgw=nbgw-ibgw+1
    nvelgw=chgval-2.d0*dble(ibgw-1)
!
!   set nbpol, the number of bands used for polmat calculations 
!
    if(emaxpol.le.0.0d0 )then 
      nbpol=nstfv
    else
      nbpol=0
      do ik=1,nkpt
        do ist=1,nstfv
          if(evaldft(ist,ik).le.emaxpol)then 
            nbpol=nbpol+1
          endif 
        enddo 
      enddo
    endif 

    write(fgw,101)'Number of bands (lapw):                        ', nstfv
    write(fgw,101)'Maximun Number of bands used in polarization:  ', nbpol
    write(fgw,101)'Number of gw bands (gw output):                ', nbandsgw
    write(fgw,*)  'Range of GW bands:                             ', ibgw,nbgw
    write(fgw,*)  'Number of valence electrons:                   ', int(chgval) 
    write(fgw,*)  'Number of valence electrons included in gw band', int(nvelgw) 
101 format(a,2i4)

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
