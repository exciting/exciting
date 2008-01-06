
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
    ! allocate contracted coefficients array
    if (allocated(apwdlm)) deallocate(apwdlm)
    allocate(apwdlm(nstsv,apwordmax,lmmaxapwtd,natmtot))
    if (allocated(apwdlm0)) deallocate(apwdlm0)
    allocate(apwdlm0(nstsv,apwordmax,lmmaxapwtd,natmtot))

    ! delete timing information of previous runs
    call filedel(trim(fnetim))

    ! write information
    write(unitout,'(a,i6)') 'Info('//thisnam//'): number of G+q vectors:', &
         ngq(iq)

    ! loop over k-points
    do ik=1,nkpt
       
       if (emattype.eq.0) then
          ! v-c and c-v band combinations
          call ematbdlims(1)
          call ematqk2(iq,ik)
       end if

       ! writing of matrix elements to file
       !+++++++++++++++++++++++++++++++++++

       ! end loop over k-points
    end do

  end subroutine ematq2

end module m_ematq2

