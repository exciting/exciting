
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: calcpmatgw
! !INTERFACE:
subroutine calcpmatgw
! !USES:
    use modinput
    use modmain, only: nkpt, ngkmax, apwordmax, lmmaxapw, natmtot, &
   &  nmatmax, nstfv, nstsv, nlotot, nlomax, lolmax, task, vkl, vgkl, &
   &  ngk, gkc, tpgkc, sfacgk, igkig, vgkc, filext
    use modxs
    use modgw
! !DESCRIPTION:
!   Calculates the momentum matrix elements using routines {\tt genpmatxs} and
!   {\tt genpmatcor}. The matrices are written then to direct access files 
!   {\tt PMAT.OUT} and {\tt PMATCOR.OUT}.
!
! !REVISION HISTORY:
!
!   Created based on writepmatxs.f90, July 2011 (DIN)
!
! !LOCAL VARIABLES:
    implicit none

    integer(4) :: ik,recl,recl2
    
    real(8)    :: tstart, tend

    complex(8), allocatable :: apwalmt(:,:,:,:)
    complex(8), allocatable :: evecfvt(:,:)
    complex(8), allocatable :: evecsvt(:,:)
    complex(8), allocatable :: pmat(:,:,:)
    complex(8), allocatable :: pmatc(:,:,:,:)

! !EXTERNAL ROUTINES: 
    external pmatrad
    external pmatradcor
    external getevecfv
    external getevecsv
    external match
    external genapwcmt
    external genlocmt
    external genpmatxs
    external genpmatxscor

!EOP
!BOC
    
    call cpu_time(tstart)
    if(tstart.lt.0.0d0)write(fgw,*)'warning, tstart < 0'
    
!   local arrays
    allocate(apwalmt(ngkmax,apwordmax,lmmaxapw,natmtot))
    allocate(evecfvt(nmatmax,nstfv))
    allocate(evecsvt(nstsv,nstsv))
      
!   allocate the momentum matrix array
    allocate(pmat(3,nstsv,nstsv))
    allocate(pmatc(3,nstsv,nclm,natmtot))

!   record length for momentum matrix elements file
    recl=16*(3*nstsv*nstsv)
    open(50,file='PMAT.OUT',action='WRITE',form='UNFORMATTED', &
   &    access='DIRECT',status='REPLACE',recl=recl)
    if(wcore)then 
      recl2=16*(3*nstsv*nclm*natmtot)
      open(51,file='PMATCOR.OUT',action='WRITE',form='UNFORMATTED', &
     &    access='DIRECT',status='REPLACE',recl=recl2)
    endif   

!----------------------------------------------------------------------!
!   calculate gradient of radial functions times spherical harmonics   !
!----------------------------------------------------------------------!

!   local arrays used for calculating gradient of radial functions
    if (allocated(apwcmt)) deallocate(apwcmt)
    allocate(apwcmt(nstsv,apwordmax,lmmaxapw,natmtot))
    if (allocated(ripaa)) deallocate(ripaa)
    allocate(ripaa(apwordmax,lmmaxapw,apwordmax,lmmaxapw,natmtot,3))
    if (nlotot .gt. 0) then
       if (allocated(locmt)) deallocate(locmt)
       allocate(locmt(nstsv,nlomax,-lolmax:lolmax,natmtot))
       if (allocated(ripalo)) deallocate(ripalo)
       allocate(ripalo(apwordmax,lmmaxapw,nlomax,-lolmax:lolmax,natmtot,3))
       if (allocated(riploa)) deallocate(riploa)
       allocate(riploa(nlomax,-lolmax:lolmax,apwordmax,lmmaxapw,natmtot,3))
       if (allocated(riplolo)) deallocate(riplolo)
       allocate(riplolo(nlomax,-lolmax:lolmax,nlomax,-lolmax:lolmax,natmtot,3))
    end if

!   valence-valence
    Call pmatrad

!   core-valence
    if (wcore) then
       if (allocated(ripacor)) deallocate(ripacor)
       allocate(ripacor(apwordmax,lmmaxapw,ncmax,nclm,natmtot,3))
       if (allocated(ripcora)) deallocate(ripcora)
       allocate(ripcora(ncmax,nclm,apwordmax,lmmaxapw,natmtot,3))
       call pmatradcor
    end if

!---------------------------------!
!  calculate the matrix elements  !
!---------------------------------!
    do ik=1,nkpt

!      get the eigenvectors and values from file
       call getevecfv(vkl(1,ik),vgkl(1,1,1,ik),evecfvt)
       call getevecsv(vkl(1,ik),evecsvt)

!      find the matching coefficients
       call match(ngk(1,ik),gkc(1,1,ik),tpgkc(1,1,1,ik), &
      &  sfacgk(1,1,1,ik),apwalmt)

!      generate APW expansion coefficients for muffin-tin
       call genapwcmt(input%groundstate%lmaxapw,ngk(1,ik),1, &
      &  nstfv,apwalmt,evecfvt,apwcmt)

!      generate local orbital expansion coefficients for muffin-tin
       if (nlotot.Gt.0) call genlocmt(ngk(1,ik),1,nstfv,evecfvt,locmt)

!      calculate the valence-valence momentum matrix elements
       call genpmatxs(ngk(1,ik),igkig(1,1,ik),vgkc(1,1,1,ik), &
      &  evecfvt,evecsvt,pmat)

       write(50,rec=ik) pmat

       if(wcore)then

!        calculate the core-valence momentum matrix elements
         call genpmatxscor(ik,pmatc)

         write(51,rec=ik) pmatc

       endif
      
    end do

    if(wcore)close(51)   
    close(50)

    write(fgw,*)
    write(fgw,'(" Info(calcpmatgw):")')
    write(fgw,'(" Momentum matrix elements written to file PMAT.OUT and PMATCOR.OUT ")')
    write(fgw,*)
    
    call cpu_time(tend)
    if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
    call write_cputime(fgw,tend-tstart, 'CALCPMATGW')

    deallocate(apwalmt,evecfvt,evecsvt)
    deallocate(pmat,pmatc)
    deallocate(apwcmt,ripaa)
    if(nlotot.gt.0)deallocate(locmt,ripalo,riploa,riplolo)
    if(wcore)deallocate(ripacor,ripcora)

end subroutine
!EOC
