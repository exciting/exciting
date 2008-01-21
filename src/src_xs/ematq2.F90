
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
    use m_putemat2
    use m_putdevalsv2
    use m_filedel
    use m_genfilname
    implicit none
    ! arguments
    integer, intent(in) :: iq
    ! local variables
    character(*), parameter :: thisnam='ematq2'
    integer :: ik
    ! filenames
    call genfilname(basename='EMAT',iq=iq,filnam=fnemat)
    call genfilname(basename='EMAT',iq=iq,procs=procs,rank=rank,&
       filnam=fnemat_t)
    call genfilname(nodotpar=.true.,basename='EMAT_TIMING',iq=iq,&
         procs=procs,rank=rank,filnam=fnetim)
    call genfilname(basename='DEVALSV',iq=iq,filnam=fndevalsv)
    call genfilname(basename='DEVALSV',iq=iq,procs=procs,rank=rank,&
       filnam=fndevalsv_t)
    ! file extension for q-point
    call genfilname(iq=iq,setfilext=.true.)
    ! calculate k+q and G+k+q related variables
    call init1xs(qvkloff(1,iq))
    ! write G+q-vectors
    call writegqpts(iq)
    ! find highest (partially) occupied and lowest (partially) unoccupied states
    call findocclims(iq,istocc0,istocc,istunocc0,istunocc,isto0,isto,istu0,istu)
    call ematbdlims(1,nst1,istlo1,isthi1,nst2,istlo2,isthi2)
    ! generate radial integrals wrt. sph. Bessel functions
    call ematrad(iq)
    ! allocate eigenvalue and eigenvector arrays
    if (allocated(evecfv)) deallocate(evecfv)
    allocate(evecfv(nmatmax,nstfv,nspnfv))
    if (allocated(evecfv0)) deallocate(evecfv0)
    allocate(evecfv0(nmatmax0,nstfv,nspnfv))
    if (allocated(evalsv0)) deallocate(evalsv0)
    allocate(evalsv0(nstsv,nkpt))
    if (allocated(occsv0)) deallocate(occsv0)
    allocate(occsv0(nstsv,nkpt))
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
    do ik=kpari,kparf
        ! set band combinations
       call ematbdlims(2*emattype,nst1,istlo1,isthi1,nst2,istlo2,isthi2)
       allocate(xiou(nst1,nst2,ngq(iq)))
       allocate(deou(nst1,nst2))
       allocate(docc12(nst1,nst2))
       call ematqk2(iq,ik)
       if (emattype.eq.0) then
          ! all band combinations
          nst3=nstsv; nst4=nstsv
          call putemat2(iq,ik,.false.,trim(fnemat_t),x1=xiou)
          call putdevalsv2(iq,ik,.false.,trim(fndevalsv_t),e1=deou,o1=docc12)
       else
          ! v-c/c-v or v-v/c-c band combinations
          allocate(xiuo(nst1,nst2,ngq(iq)))
          allocate(deuo(nst1,nst2))
          allocate(docc21(nst1,nst2))
          xiuo(:,:,:)=xiou(:,:,:)
          deuo(:,:)=deou(:,:)
          docc21(:,:)=docc12(:,:)
          deallocate(xiou)
          deallocate(deou)
          deallocate(docc12)
          call ematbdlims(2*emattype-1,nst1,istlo1,isthi1,nst2,istlo2,isthi2)
          nst3=nst2; nst4=nst1
          allocate(xiou(nst1,nst2,ngq(iq)))
          allocate(deou(nst1,nst2))
          allocate(docc12(nst1,nst2))
          call ematqk2(iq,ik)
          call putemat2(iq,ik,.false.,trim(fnemat_t),x1=xiou,x2=xiuo)
          call putdevalsv2(iq,ik,.false.,trim(fndevalsv_t),e1=deou,o1=docc12, &
               e2=deuo,o2=docc21)
       end if
       deallocate(xiou)
       deallocate(deou)
       deallocate(docc12)
       if (allocated(xiuo)) deallocate(xiuo)
       if (allocated(deuo)) deallocate(deuo)
       if (allocated(docc21)) deallocate(docc21)
       ! synchronize for common number of k-points to all processes
       if ((partype.eq.'k').and.(ik-kpari+1 <= nkpt/procs)) call barrier
       ! end loop over k-points
    end do
  end subroutine ematq2

end module m_ematq2

