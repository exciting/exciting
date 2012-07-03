
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
    integer(4)              :: recl,recl2,ik
    complex(8), allocatable :: apwalmt(:,:,:,:)
    complex(8), allocatable :: evecfvt(:,:)
    complex(8), allocatable :: evecsvt(:,:)
    complex(8), allocatable :: pmat(:,:,:)
    complex(8), allocatable :: pmatc(:,:,:)
    
    real(8) :: tstart, tend
    
    call cpu_time(tstart)
    if(tstart.lt.0.0d0)write(fgw,*)'warning, tstart < 0'
    
    allocate(apwalmt(ngkmax,apwordmax,lmmaxapw,natmtot))
    allocate(evecfvt(nmatmax,nstfv))
    allocate(evecsvt(nstsv,nstsv))
      
!   allocate the momentum matrix array
    allocate(pmat(3,nstfv,nstfv))
    if(wcore)then
      allocate(pmatc(3,ncg,nstfv))
    end if

!   record length for momentum matrix elements file
    recl=16*(3*nstsv*nstsv)
    open(50,file='PMAT.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
     status='REPLACE',recl=recl)
    if(wcore)then 
      recl2=16*(3*ncg*nstsv)
      open(51,file='PMATCOR.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
       status='REPLACE',recl=recl2)
    endif   

    do ik=1,nkpt

!     get the eigenvectors and values from file
      call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfvt)
      call getevecsv(vkl(:,ik),evecsvt)
      
!     find the matching coefficients
      call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik), &
     &    sfacgk(:,:,1,ik),apwalmt)

!     valence-valence matrix elements
      call genpmat(ngk(1,ik),igkig(1,1,ik),vgkc(1,1,1,ik),apwalmt, &
     &  evecfvt,evecsvt,pmat)
      write(50,rec=ik) pmat

!     core-valence contribution      
      if(wcore)then
        call genpmatcor(ik,ngk(1,ik),apwalmt,evecfvt,evecsvt,pmatc)
        write(51,rec=ik) pmatc
      endif

    end do

    if(wcore)close(51)   
    close(50)

    write(fgw,*)
    write(fgw,'(" Info(calcpmat):")')
    write(fgw,'(" Momentum matrix elements written to file PMAT.OUT and PMATCOR.OUT ")')
    write(fgw,*)
    
    call cpu_time(tend)
    if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
    call write_cputime(fgw,tend-tstart, 'CALCPMAT')

    deallocate(apwalmt,evecfvt,evecsvt,pmat)
    if(wcore)then
      deallocate(pmatc)
    end if
end subroutine
!EOC
