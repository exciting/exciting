
! Copyright (C) 2006-2007 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_ematq2
  implicit none
contains

  subroutine ematq2(iq)
    use modmain
    use modxs
    use modmpi
    use m_ematrad
    use m_ematqk2
    use m_writegqpts
    use m_getunit
    use m_filedel
    use m_genfilname
    implicit none
    ! arguments
    integer, intent(in) :: iq
    ! local variables
    character(*), parameter :: thisnam = 'ematq2'
    integer :: ik,ki,kf

    ! filenames
    call genfilname(basename='EMAT',etype=emattype,iq=iq,filnam=fnemat)
    call genfilname(basename='EMAT',etype=emattype,iq=iq,procs=procs,rank=rank,&
       filnam=fnemat_t)
    call genfilname(nodotpar=.true.,basename='EMAT_TIMING',iq=iq,&
         procs=procs,rank=rank,filnam=fnetim)

    ! set limits for band combinations
    call ematbdlims(1)

    ! write G+q-vectors
    call writegqpts(iq)

    ! generate radial integrals wrt. sph. Bessel functions
    call ematrad(iq)

    ! allocate eigenvalue and eigenvector arrays
    if (allocated(evecfv)) deallocate(evecfv)
    if (allocated(evecfv0)) deallocate(evecfv0)
    allocate(evecfv(nmatmax,nstfv,nspnfv))
    allocate(evecfv0(nmatmax0,nstfv,nspnfv))

    ! allocate helper matrix
    if (allocated(xih)) deallocate(xih)
    allocate(xih(nlotot,nlotot))
    ! allocate matrix elements array
    if (allocated(xiou)) deallocate(xiou)
    allocate(xiou(nst1,nst2,ngq(iq)))

    if (allocated(xiohalo)) deallocate(xiohalo)
    if (allocated(xiuhloa)) deallocate(xiuhloa)
    if (allocated(xihlolo)) deallocate(xihlolo)
    allocate(xiohalo(nst1,nlotot))
    allocate(xiuhloa(nlotot,nst2))
    allocate(xihlolo(nlotot,nlotot))
    if (allocated(apwdlm)) deallocate(apwdlm)
    if (allocated(apwdlm0)) deallocate(apwdlm0)
    allocate(apwdlm(nstsv,apwordmax,lmmaxapwtd,natmtot))
    allocate(apwdlm0(nstsv,apwordmax,lmmaxapwtd,natmtot))

    ! delete timing information of previous runs
    call filedel(trim(fnetim))

    ! write information
    write(unitout,'(a,i6)') 'Info('//thisnam//'): number of G+q vectors:', &
         ngq(iq)

    ! limits for k-point loop
    ki=kpari
    kf=kparf

!!$    call getunit(un)
    ! loop over k-points
    do ik = ki, kf
       call ematqk2(iq,ik)
!!$#ifdef MPI
!!$       if (ik-ki+1 <= nkpt/procs) then
!!$          ! synchronize for common number of k-points to all processes
!!$          call barrier(rank=rank,procs=procs,un=un,async=0,string='.barrier')
!!$       end if
!!$#endif
       ! end loop over k-points
    end do

  end subroutine ematq2

end module m_ematq2

