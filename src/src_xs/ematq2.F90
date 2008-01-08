
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
    use m_genfilname
    use m_writegqpts
    use m_ematrad
    use m_filedel
    use m_ematqk2
    use m_putemat2
    implicit none
    ! arguments
    integer, intent(in) :: iq
    ! local variables
    character(*), parameter :: thisnam = 'ematq2'
    integer :: ik
    ! filenames
    call genfilname(basename='EMAT',etype=emattype,iq=iq,filnam=fnemat)
    call genfilname(basename='EMAT_TIMING',etype=emattype,iq=iq,filnam=fnetim)
    ! write G+q-vectors
    call writegqpts(iq)
    ! generate radial integrals wrt. sph. Bessel functions
    call ematrad(iq)
    ! allocate eigenvalue and eigenvector arrays
    if (allocated(evecfv)) deallocate(evecfv)
    if (allocated(evecfv0)) deallocate(evecfv0)
    allocate(evecfv(nmatmax,nstfv,nspnfv))
    allocate(evecfv0(nmatmax,nstfv,nspnfv))
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
       ! set band combinations
       call ematbdlims(2*emattype,nst1,istlo1,isthi1,nst2,istlo2,isthi2)
       allocate(xiou(nst1,nst2,ngq(iq)))
       call ematqk2(iq,ik)
       if (emattype.eq.0) then
          ! all band combinations
          call putemat2(iq,ik,trim(fnemat),xiou)
       else
          ! v-c/c-v or v-v/c-c band combinations
          allocate(xiuo(nst1,nst2,ngq(iq)))
          xiuo(:,:,:)=xiou(:,:,:)
          deallocate(xiou)
          call ematbdlims(2*emattype-1,nst1,istlo1,isthi1,nst2,istlo2,isthi2)
          allocate(xiou(nst1,nst2,ngq(iq)))
          call ematqk2(iq,ik)
          call putemat2(iq,ik,trim(fnemat),xiou,xiuo)          
       end if
       

!!$       select case(emattype)
!!$       case(0)
!!$          ! all band combinations
!!$          call ematbdlims(0,nst1,istlo1,isthi1,nst2,istlo2,isthi2)
!!$          allocate(xiou(nst1,nst2,ngq(iq)))
!!$          call ematqk2(iq,ik)
!!$          ! write to file
!!$          call putemat2(iq,ik,trim(fnemat),xiou)
!!$       case(1)
!!$          ! c-v band combinations
!!$          call ematbdlims(2,nst1,istlo1,isthi1,nst2,istlo2,isthi2)
!!$          allocate(xiou(nst1,nst2,ngq(iq)))
!!$          call ematqk2(iq,ik)
!!$          ! save to array xiuo
!!$          allocate(xiuo(nst1,nst2,ngq(iq)))
!!$          xiuo(:,:,:)=xiou(:,:,:)
!!$          deallocate(xiou)
!!$          ! v-c band combinations
!!$          call ematbdlims(1,nst1,istlo1,isthi1,nst2,istlo2,isthi2)
!!$          allocate(xiou(nst1,nst2,ngq(iq)))
!!$          call ematqk2(iq,ik)
!!$          ! write to file
!!$          call putemat2(iq,ik,trim(fnemat),xiou,xiuo)
!!$       case(2)
!!$          ! v-v band combinations
!!$          call ematbdlims(4,nst1,istlo1,isthi1,nst2,istlo2,isthi2)
!!$          allocate(xiou(nst1,nst2,ngq(iq)))
!!$          call ematqk2(iq,ik)
!!$          ! save to array xiuo
!!$          allocate(xiuo(nst1,nst2,ngq(iq)))
!!$          xiuo(:,:,:)=xiou(:,:,:)
!!$          deallocate(xiou)
!!$          ! c-c band combinations
!!$          call ematbdlims(3,nst1,istlo1,isthi1,nst2,istlo2,isthi2)
!!$          allocate(xiou(nst1,nst2,ngq(iq)))
!!$          call ematqk2(iq,ik)
!!$          ! write to file
!!$          call putemat2(iq,ik,trim(fnemat),xiou,xiuo)
!!$       end select


       deallocate(xiou)
       if (allocated(xiuo)) deallocate(xiuo)
       ! end loop over k-points
    end do
  end subroutine ematq2

end module m_ematq2

