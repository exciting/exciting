!BOP
!
! !ROUTINE: kintw
!
! !INTERFACE:
      subroutine kintw()
      
! !DESCRIPTION:
!
! This subroutine calculates the weights for BZ integrations.
!
!
! !USES:
!
      use modmain
      use modgw
     
! !LOCAL VARIABLES:

      implicit none

      integer(4) :: ia
      integer(4) :: ias
      integer(4) :: is
      integer(4) :: ist
      integer(4) :: ik, ikp
      logical    :: lprt=.false.
      real(8)    :: tstart, tend
      
      real(8), allocatable :: bandpar(:,:) ! Band for which the
!                                            weights are being calculated
      real(8), allocatable :: cwpar(:,:)   ! weight of one band
!
!EOP
!BOC
      call cpu_time(tstart)
      if(tstart.lt.0.0d0)write(98,*)'warning, tstart < 0'

      if (allocated(kiw)) deallocate(kiw)
      allocate(kiw(nstfv,nkptnr))
      if (allocated(kwfer)) deallocate(kwfer)
      allocate(kwfer(nstfv,nkptnr))
      kiw=0.0d0
      kwfer=0.0d0
      
      if (iopcore<2) then
        if (allocated(ciw)) deallocate(ciw)
        allocate(ciw(natmtot,ncmax))
        ciw=0.0d0
      end if
!
!     Q-dependent integration weights used to calculate the polarizability
!     (qdepw.f90)
!
      if (allocated(kcw)) deallocate(kcw)
      allocate(kcw(nstfv,nstfv,1:nomeg,nkptnr))
      if (iopcore<2) then
        if (allocated(unw)) deallocate(unw)
        allocate(unw(natmtot,ncmax,nstfv,1:nomeg,nkptnr))
      end if
      if((testid.eq.7).and.(fflg.eq.2)) then
        if (allocated(kcwsurf)) deallocate(kcwsurf)
        allocate(kcwsurf(nstfv,nstfv,1:nomeg,nkptnr))
        kcwsurf=0.0d0
        if (iopcore<2) then
          if (allocated(unwsurf)) deallocate(unwsurf)
          allocate(unwsurf(natmtot,ncmax,nstfv,1:nomeg,nkptnr))
          unwsurf=0.0d0
        end if
      endif
!
!     allocate local arrays
!
      metallic = .false.
      
      if (iopcore<2) then
!---------------------------------------------------------------------
!                 core
!---------------------------------------------------------------------
        allocate(bandpar(1,nkptnr))
        allocate(cwpar(1,nkptnr))
        do is=1,nspecies
          do ia=1,natoms(is)
            ias=idxas(ia,is)
!           Product of core and apw functions
            do ist=1,ncore(is)
              bandpar(1,:)=evalcr(ist,ias)       
              call tetiw(nkptnr,ntetnr,1,bandpar,tnodesnr,wtetnr,tvol,efermi,cwpar)
              ciw(ias,ist)=cwpar(1,1)
            enddo ! ist
          enddo ! ia
        enddo ! is
        deallocate(bandpar,cwpar)
      end if ! iopcore

!---------------------------------------------------------------------
!                 valence
!---------------------------------------------------------------------
      allocate(bandpar(nstfv,nkptnr))
      do ik=1,nkptnr
         ikp=indkp(ik)
         bandpar(1:nstfv,ik)=evaldft(1:nstfv,ikp)
      enddo  

      call tetiw(nkptnr,ntetnr,nstfv,bandpar,tnodesnr,wtetnr,tvol,efermi,kiw)

      call tetiwsurf(nkptnr,ntetnr,nstfv,bandpar,tnodesnr,wtetnr,tvol,efermi,kwfer)
      
      deallocate(bandpar)
      
      if (maxval(kwfer).ne.0.0d0) metallic=.true.
      !write(fgw,*)'kwfer: metallic =', metallic
      
      if(lprt)then
        if (iopcore<2) then
          write(15,*)'core: ciw'
          write(15,*) ciw(1,:)
        end if
        write(15,*)'valence: kiw'
        do ik=1,nkptnr
        write(15,*) ik, kiw(:,ik)
        enddo
        stop
      endif
      
      call cpu_time(tend)
      if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
      call write_cputime(fgw,tend-tstart, 'KINTW')
      
      return
      end subroutine kintw
!EOC          
         
