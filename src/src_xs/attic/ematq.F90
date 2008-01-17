
! Copyright (C) 2006-2007 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_ematq
  implicit none
contains

  subroutine ematq(iq)
    use modmain
    use modxs
    use modmpi
    use m_ematrad
    use m_ematqk
    use m_writegqpts
    use m_filedel
    use m_genfilname
    implicit none
    ! arguments
    integer, intent(in) :: iq
    ! local variables
    character(*), parameter :: thisnam = 'ematq'
    integer :: ik
    real(8) :: vkloff_save(3)

    ! filenames
    call genfilname(basename='EMAT',iq=iq,filnam=fnemat)
    call genfilname(basename='EMAT',iq=iq,procs=procs,rank=rank,&
       filnam=fnemat_t)
    call genfilname(nodotpar=.true.,basename='EMAT_TIMING',iq=iq,&
         procs=procs,rank=rank,filnam=fnetim)
    call genfilname(basename='DEVALSV',iq=iq,filnam=fndevalsv)
    call genfilname(basename='DEVALSV',iq=iq,procs=procs,rank=rank,&
       filnam=fndevalsv_t)

    ! save k-point offset
    vkloff_save(:)=vkloff(:)

    ! file extension for q-point
    call genfilname(iq=iq,setfilext=.true.)
    ! shift k-mesh by q-point
    vkloff(:)=qvkloff(:,iq)

    ! calculate k+q and G+k+q related variables
    call init1xs

    ! write G+q-vectors
    call writegqpts(iq)

    ! find highest (partially) occupied and lowest (partially) unoccupied states
    call findocclims(iq,istocc0,istocc,istunocc0,istunocc,isto0,isto,istu0,istu)
    call ematbdlims(1,nst1,istlo1,isthi1,nst2,istlo2,isthi2)

    ! generate radial integrals wrt. sph. Bessel functions
    call ematrad(iq)

    ! allocate eigenvalue and eigenvector arrays
    if (allocated(evecfv)) deallocate(evecfv)
    if (allocated(evecfv0)) deallocate(evecfv0)
    if (allocated(evalsv0)) deallocate(evalsv0)
    allocate(evecfv(nmatmax,nstfv,nspnfv))
    allocate(evecfv0(nmatmax0,nstfv,nspnfv))
    allocate(evalsv0(nstsv,nkpt))

    ! allocate arrays for eigenvalue differences
    if (allocated(deou)) deallocate(deou)
    if (allocated(deuo)) deallocate(deuo)
    allocate(deou(nstval,nstcon))
    allocate(deuo(nstcon,nstval))
    if (allocated(docc12)) deallocate(docc12)
    allocate(docc12(nst1,nst2))
    if (allocated(docc21)) deallocate(docc21)
    allocate(docc21(nst2,nst1))
    if (allocated(occsv0)) deallocate(occsv0)
    allocate(occsv0(nstsv,nkpt))
    ! allocate helper matrix
    if (allocated(xih)) deallocate(xih)
    allocate(xih(nlotot,nlotot))
    ! allocate matrix elements array
    if (allocated(xiou)) deallocate(xiou)
    if (allocated(xiuo)) deallocate(xiuo)
    allocate(xiou(nstval,nstcon,ngq(iq)))
    allocate(xiuo(nstcon,nstval,ngq(iq)))

    if (allocated(xiohalo)) deallocate(xiohalo)
    if (allocated(xiuhalo)) deallocate(xiuhalo)
    if (allocated(xiohloa)) deallocate(xiohloa)
    if (allocated(xiuhloa)) deallocate(xiuhloa)
    if (allocated(xihlolo)) deallocate(xihlolo)
    allocate(xiohalo(nstval,nlotot))
    allocate(xiuhalo(nstcon,nlotot))
    allocate(xiohloa(nlotot,nstval))
    allocate(xiuhloa(nlotot,nstcon))
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

    ! loop over k-points
    do ik=kpari,kparf
       call ematqk(iq,ik)


!?????????????????????????????????????????????????????????????????????



       ! synchronize for common number of k-points to all processes
       if (ik-kpari+1 <= nkpt/procs) call barrier
       ! end loop over k-points
    end do

    ! restore offset
    vkloff(:)=vkloff_save(:)
    ! restore file extension
    call genfilname(setfilext=.true.)

  end subroutine ematq

end module m_ematq

