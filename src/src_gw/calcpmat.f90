
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: calcpmat
! !INTERFACE:
subroutine calcpmat
! !USES:
    use modmain
    use modgw
    use modxs
    use modmpi
    use m_getunit
! !DESCRIPTION:
!   Calculates the momentum matrix elements using routine {\tt genpmat} and
!   writes them to direct access file {\tt PMAT.OUT}.
!
! !REVISION HISTORY:
!   Created from writepmat August 2006 (RGA)
!   Revisited: June 2011 by DIN
!
!EOP
!BOC
    implicit none
    integer(4)              :: ik,ik0
    integer(4)              :: fid, Recl
    complex(8), allocatable :: apwalmt(:,:,:,:)
    complex(8), allocatable :: evecfvt(:,:)
    complex(8), allocatable :: evecsvt(:,:)
    complex(8), allocatable :: pmat(:,:,:,:)
    complex(8), allocatable :: pmatc(:,:,:,:)
    real(8) :: tstart, tend
    integer :: ikstart, ikend
    integer, allocatable :: ikp2rank(:)
    
    character(128)::sbuffer
    call cpu_time(tstart)
    if(tstart.lt.0.0d0)write(fgw,*)'warning, tstart < 0'
    
#ifdef MPI    
    ikstart = firstofset(rank,nkpt)
    ikend = lastofset(rank,nkpt)
#else
    ikstart = 1
    ikend = nkpt
#endif

    allocate(ikp2rank(nkpt))
    ikp2rank = -1
    do ik = ikstart, ikend
      ikp2rank(ik) = rank
    end do ! ik
    
    allocate(apwalmt(ngkmax,apwordmax,lmmaxapw,natmtot))
    allocate(evecfvt(nmatmax,nstfv))
    allocate(evecsvt(nstsv,nstsv))
      
!   allocate the momentum matrix array

    allocate(pmat(3,nstfv,nstfv,ikstart:ikend))
    pmat = zzero
    if (iopcore.eq.0) then 
      allocate(pmatc(3,ncg,nstfv,ikstart:ikend))
      pmatc = zzero
    endif   

    do ik = ikstart, ikend

      ik0=idikp(ik)

!     get the eigenvectors and values from file
      call getevecfvgw(ik0,evecfvt)
      call getevecsvgw(ik0,evecsvt)
      
!     find the matching coefficients
      call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik), &
     &  sfacgk(:,:,1,ik),apwalmt)

!     valence-valence matrix elements
      call genpmat(ngk(1,ik),igkig(1,1,ik),vgkcnr(1,1,1,ik), &
     &  apwalmt,evecfvt,evecsvt,pmat(:,:,:,ik))
 
!     core-valence contribution      
      if (iopcore.eq.0) then
        call genpmatcor(ik,ngk(1,ik),apwalmt,evecfvt,evecsvt,pmatc(:,:,:,ik))
      endif

    end do
    
    deallocate(apwalmt,evecfvt,evecsvt)
    
    !==========================
    ! Write results to files
    !==========================
    
    ! overwrite existing files
    if (rank==0) then
      call getunit(fid)
      open(fid,File='PMAT.OUT',form='UNFORMATTED',status='REPLACE')
      close(fid)
      if (input%gw%coreflag=='all') then
        call getunit(fid)
        open(fid,File='PMATCOR.OUT',form='UNFORMATTED',status='REPLACE')
        close(fid)
      end if
    endif
    call barrier
    
    do ik = 1, nkpt
      if (rank==ikp2rank(ik)) then
        call getunit(fid)
        inquire(iolength=recl) pmat(:,:,:,ik)
        open(fid,File='PMAT.OUT',Action='WRITE',Form='UNFORMATTED',&
        &    Access='DIRECT',Status='OLD',Recl=recl)
        write(fid,rec=ik) pmat(:,:,:,ik)
        close(fid)
        if (iopcore.eq.0) then
          call getunit(fid)
          inquire(iolength=recl) pmatc(:,:,:,ik)
          open(fid,File='PMATCOR.OUT',Action='WRITE',Form='UNFORMATTED',&
          &    Access='DIRECT',Status='OLD',Recl=recl)
          write(fid,rec=ik) pmatc(:,:,:,ik)
          close(fid)
        end if
      end if ! rank
      call barrier
    end do ! ikp

    deallocate(pmat)
    if (iopcore.eq.0) deallocate(pmatc)

    if (rank==0) then
      write(fgw,*)
      write(fgw,'(" Info(calcpmat):")')
      write(fgw,'(" Momentum matrix elements written to file PMAT.OUT and PMATCOR.OUT ")')
      write(fgw,*)
    end if
    
    call cpu_time(tend)
    if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
    call write_cputime(fgw,tend-tstart, 'CALCPMAT')

end subroutine
!EOC
