
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
    use m_getunit
    use m_filedel
    use m_genfilname
    implicit none
    ! arguments
    integer, intent(in) :: iq
    ! local variables
    character(*), parameter :: thisnam = 'ematq'
    integer :: ik,un,ki,kf
    real(8) :: vkloff_save(3)


!!$integer :: iv,ic,isym,lspl,igq,iklt(3),ikt,igqt,ivl(3),nsym
!!$real(8) :: frtc(3),c(3,3),vl(3),vc(3),vct(3),kct(3),klt(3),rklt(3),pklt,t1
!!$real(8) :: vg(3),vgcc(3),vgct(3),vglt(3)
!!$complex(8) :: zt1,zt2


    ! filenames
    call genfilname(basename='EMAT',iq=iq,filnam=fnemat)
    call genfilname(basename='EMAT',iq=iq,procs=procs,rank=rank,&
       filnam=fnemat_t)
    call genfilname(nodotpar=.true.,basename='EMAT_TIMING',iq=iq,&
         procs=procs,rank=rank,filnam=fnetim)
    call genfilname(basename='DEVALSV',iq=iq,filnam=fndevalsv)
    call genfilname(basename='DEVALSV',iq=iq,procs=procs,rank=rank,&
       filnam=fndevalsv_t)

    ! checkpointing starts
    call getunit(un)
!!$    if (.not.tresume) then
!!$       ! k-point index
!!$       resumechkpts(1,1)=0
!!$       call resupd(un,task,resumechkpts,' : k-point index')
!!$    else
!!$       ! jump into next checkpoint (k-point)
!!$       resumechkpts(1,1)=resumechkpts(1,1)+1
!!$    end if

    ! save k-point offset
    vkloff_save = vkloff

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
!!$    if (tresume) ki=resumechkpts(1,1)
    call getunit(un)
    ! loop over k-points
    do ik = ki, kf
       call ematqk(iq,ik)




!!$nsym=nsymcrysstr(ik)
!!$do iv=1,nstval
!!$   do ic=1,nstcon
!!$      do igq=1,ngq(iq)
!!$         vl(:)=vgql(:,igq,iq)
!!$         vc=matmul(bvec,vl)
!!$         do isym=1,nsymcrys
!!$            lspl=lsplsymc(isym)
!!$            c(:,:)=symlatc(:,:,lspl)
!!$            vct=matmul(c,vc)
!!$            frtc(:)=matmul(avec,vtlsymc(:,isym))
!!$            t1=dot_product(vct,frtc)
!!$            zt1=cmplx(cos(t1),sin(t1),8)
!!$
!!$            kct=matmul(c,vkc(:,ik))
!!$            klt=matmul(binv,kct)
!!$            call r3frac(epslat,klt,ivl)
!!$            rklt=(klt-vkloff/ngridk)*ngridk
!!$            iklt=nint(rklt)
!!$            ikt=ikmapnr(iklt(1),iklt(2),iklt(3))
!!$            pklt=sum(abs(rklt-iklt)*ngridk)*100.d0
!!$
!!$            vg(:)=dble(ivg(:,igqig(igq,iq)))
!!$            vgcc=matmul(bvec,vg)
!!$            vgct=matmul(c,vgcc)
!!$            vglt=matmul(binv,vgct)
!!$            ivl=nint(vglt)
!!$
!!$            igqt=ivgigq(ivl(1),ivl(2),ivl(3),iq)
!!$            zt2=xiou(iv,ic,igqt)
!!$
!!$            write(1600,'(5i4,2f8.2,3x,2f8.2,2i4)') &
!!$                 ik,iv,ic,igq,isym,&
!!$                 xiou(iv,ic,igq),zt1*zt2,igqt,ikt
!!$         end do
!!$      end do
!!$   end do
!!$end do






!!$       resumechkpts(1,1)=ik
!!$       call resupd(un,task,resumechkpts,' : k-point index')
#ifdef MPI
       if (ik-ki+1 <= nkpt/procs) then
          ! synchronize for common number of k-points to all processes
          call barrier(rank=rank,procs=procs,un=un,async=0,string='.barrier')
       end if
#endif
       ! end loop over k-points
    end do

    ! close files
    close(unit1)

    ! restore offset
    vkloff = vkloff_save
    ! restore file extension
    call genfilname(setfilext=.true.)

!!$    if (tresume) tresume=.false.

  end subroutine ematq

end module m_ematq

