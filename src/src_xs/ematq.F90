
! Copyright (C) 2006-2007 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_ematq
  implicit none
contains

  subroutine ematq(iq)
    use modmain
    use modtddft
    use modpar
    use m_ematrad
    use m_ematqk
    use m_writegqpts
    use m_getunit
    use m_filedel
    use m_genfilname
    implicit none
    ! arguments
    integer, intent(in) :: iq
    ! local variables
    character(*), parameter :: thisnam = 'ematq'
    integer :: ik,un,recl,ki,kf
    real(8) :: stim, vkloff_save(3)

    ! filenames
    call genfilname(basename='EMAT',iq=iq,filnam=fnemat)
    call genfilname(basename='EMAT',iq=iq,nproc=nproc,rank=rank-1,&
       filnam=fnemat_t)
    call genfilname(nodotpar=.true.,basename='EMAT_TIMING',iq=iq,&
         nproc=nproc,rank=rank-1,filnam=fnetim)
    call genfilname(basename='DEVALSV',iq=iq,filnam=fndevalsv)
    call genfilname(basename='DEVALSV',iq=iq,nproc=nproc,rank=rank-1,&
       filnam=fndevalsv_t)

    ! checkpointing starts
    call getunit(un)
    if (.not.tresume) then
       ! k-point index
       resumechkpts(1,1)=0
       call resupd(un,task,resumechkpts,' : k-point index')
    else
       ! jump into next checkpoint (k-point)
       resumechkpts(1,1)=resumechkpts(1,1)+1
    end if

    ! save k-point offset
    vkloff_save = vkloff

    ! file extension for q-point
    call genfilname(iq=iq,setfilext=.true.)
    ! shift k-mesh by q-point
    vkloff(:)=qvkloff(:,iq)
    ! calculate k+q and G+k+q related variables
    call init1td

    ! write G+q-vectors
    call writegqpts(iq)

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
    if(allocated(deou)) deallocate(deou)
    if(allocated(deuo)) deallocate(deuo)
    allocate(deou(nstval,nstcon))
    allocate(deuo(nstcon,nstval))
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

    ! limits for k-point loop
    ki=kpari
    kf=kparf
    ! resume task, first checkpoint is k-point index
    if (tresume) ki=resumechkpts(1,1)
    call getunit(un)
    ! loop over k-points
    do ik = ki, kf
       call ematqk(iq,ik)
       resumechkpts(1,1)=ik
       call resupd(un,task,resumechkpts,' : k-point index')
#ifdef MPI
       if (ik-ki+1 <= nkpt/nproc) then
          ! synchronize for common number of k-points to all processes
          call barrier(rank=rank,nproc=nproc,un=un,async=0,string='.barrier')
       end if
#endif
    end do

    ! close files
    close(unit1)

    ! restore offset
    vkloff = vkloff_save
    ! restore file extension
    call genfilname(setfilext=.true.)

    if (tresume) tresume=.false.

  end subroutine ematq

end module m_ematq

