
! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: dfq
! !INTERFACE:
subroutine dfq(iq)
! !USES:
  use modmain
  use modxs
  use modtetra
  use modmpi
  use m_genwgrid
  use m_getpemat
  use m_dftim
  use m_gettetcw
  use m_putx0
  use m_getunit
  use m_writevars
  use m_filedel
  use m_genfilname
! !DESCRIPTION:
!   Calculates the symmetrized Kohn-Sham response function $\chi^0_{\bf{GG'}}
!   ({\bf q},\omega)$ for one ${\bf q}$-point according to
!   $$  \chi^0_{\bf{GG'}}({\bf q},\omega) = \sum_{ou{\bf k}} \left[
!      M^{\bf G}_{ou{\bf k}}({\bf q}) M^{\bf G'}_{ou{\bf k}}({\bf q})^*
!      w_{ou{\bf k}}({\bf q},\omega) +
!      M^{\bf G}_{uo{\bf k}}({\bf q}) M^{\bf G'}_{uo{\bf k}}({\bf q})^*
!      w_{uo{\bf k}}({\bf q},\omega) \right]
!   $$
!   It is related to the Kohn-Sham response function
!    $\chi^0_{\bf{GG'}}({\bf q},\omega)$ by
!   $$ \chi^0_{\bf{GG'}}({\bf q},\omega) = v^{-\frac{1}{2}}_{\bf G}({\bf q})
!      \bar{\chi}^0_{\bf{GG'}}({\bf q},\omega)
!      v^{-\frac{1}{2}}_{\bf G'}({\bf q}) $$
!   and is well defined in the limit as ${\bf q}$ goes to zero.
!   The symmetrized matrix elements are defined as
!   $$   M^{\bf G}_{ou{\bf k}}({\bf q}) =
!         v^{-\frac{1}{2}}_{\bf G}({\bf q})
!         \bar{M}^{\bf G}_{ou{\bf k}}({\bf q}), $$
!   where
!   $$   \bar{M}^{\bf G}_{ou{\bf k}}({\bf q}) =
!        \langle o{\bf k}|e^{-i({\bf{G+q}}){\bf r}}|u{\bf k+q} \rangle. $$
!   For ${\bf G}=0$ we have to consider three vectors stemming from the limits
!   as ${\bf q}\rightarrow 0$ along the Cartezian basis vectors ${\bf e_i}$,
!   i.e., we can think of ${\bf 0}_1,{\bf 0}_2,{\bf 0}_3$ in place of ${\bf 0}$.
!   The weights $w_{\rm ou{\bf k}}({\bf q},\omega)$ are defined as
!   $$   w_{nm{\bf k}}({\bf q},\omega) = \lambda_{\bf k}
!        \frac{f_{n{\bf k}}-f_{m{\bf k+q}}}
!        {\varepsilon_{n{\bf k}}-\varepsilon_{m{\bf k+q}}+
!        \Delta_{n{\bf k}}-\Delta_{m{\bf k+q}}+\omega+i\eta} $$
!   in the case where we use a Lorentzian broadening $\eta$.
!   In the above expression $\lambda_{\bf k}$ is the weight of the
!   ${\bf k}$-point, $\varepsilon_{n{\bf k}}$ and
!   $\varepsilon_{m{\bf k+q}}$ are the DFT Kohn-Sham energies,
!   $\Delta_{n{\bf k}}$ and $\Delta_{m{\bf k+q}}$ are the scissors corrections
!   that are non-zero in the case where $m{\bf k+q}$ corresponds to a
!   conduction state.
!   The indices $o$ and $u$ denote {\it at least partially occupied} and
!   {\it at least partially unoccupied} states, respectively.
!   The symmetrized Kohn-Sham response function can also be calculated
!   for imaginary frequencies $i\omega$ without broadening $\eta$. In this
!   case the replacement
!   $$ \omega+i\eta \mapsto i\omega $$
!   is applied to the expressions for the weights.
!   Optionally, the weights can be calculated with the help of the linear
!   tetrahedron method (including Bloechl's correction).
!   This routine can be run with MPI parallelization for energies.
!
! !REVISION HISTORY:
!   Created March 2005 (Sagmeister)
!   Added band and k-point analysis, 2007-2008 (Sagmeister)
!EOP
!BOC
  implicit none
  ! arguments
  integer, intent(in) :: iq
  ! local variables
  character(*), parameter :: thisnam='dfq'
  character(256) :: fnscreen
  real(8), parameter :: epstetra=1.d-8
  complex(8), allocatable :: w(:)
  complex(8), allocatable :: chi0(:,:,:),hdg(:,:,:)
  complex(8), allocatable :: chi0w(:,:,:,:),chi0h(:,:,:)
  complex(8), allocatable :: wou(:),wuo(:),wouw(:),wuow(:),wouh(:),wuoh(:)
  complex(8), allocatable :: zvou(:),zvuo(:), chi0hs(:,:,:),bsedg(:,:)
  real(8), allocatable :: wreal(:),cw(:),cwa(:),cwsurf(:)
  real(8), allocatable :: cwt(:,:),cw1k(:,:,:),cwa1k(:,:,:),cwsurf1k(:,:,:)
  real(8), allocatable :: scis12(:,:),scis21(:,:)
  real(8) :: brd,cpu0,cpu1,cpuread,cpuosc,cpuupd,cputot,r1
  integer :: n,i,j,i1,i2,j1,j2,ik,ikq,igq,iw,wi,wf,ist1,ist2,nwdfp
  integer :: oct1,oct2,un,ig1,ig2
  logical :: tq0
  logical, external :: tqgamma,transik,transijst
  if (acont.and.tscreen) then
     write(*,*)
     write(*,'("Error(",a,"): analytic continuation does not work for &
          &screening")')
     write(*,*)
     call terminate
  end if
  ! sampling of Brillouin zone
  bzsampl=0
  if (tetradf) bzsampl=1
  ! initial and final w-point
  wi=wpari
  wf=wparf
  nwdfp=wf-wi+1
  ! matrix size for response function
  n=ngq(iq)
  ! zero broadening for analytic contiunation
  brd=broad
  if (acont) brd=0.d0
  ! zero broadening for dielectric matrix (w=0) for band-gap systems
  if (task.eq.430) brd=0.d0
  ! file extension for q-point
  if (.not.tscreen) call genfilname(iqmt=iq,setfilext=.true.)
  ! filenames for output
  if (tscreen) then
     call genfilname(basename='TETW',iq=iq,appfilext=.true.,filnam=fnwtet)
     call genfilname(basename='PMAT',appfilext=.true.,filnam=fnpmat)
     call genfilname(basename='SCREEN',bzsampl=bzsampl,&
          iq=iq,filnam=fnscreen)
     call genfilname(nodotpar=.true.,basename='EMAT_TIMING',iq=iq,&
          etype=emattype,procs=procs,rank=rank,appfilext=.true.,filnam=fnetim)
     call genfilname(nodotpar=.true.,basename='X0_TIMING',iq=iq,&
          bzsampl=bzsampl,acont=acont,procs=procs,rank=rank, &
          appfilext=.true.,filnam=fnxtim)
  else
     call genfilname(basename='TETW',iqmt=iq,filnam=fnwtet)
     call genfilname(basename='PMAT_XS',filnam=fnpmat)
     call genfilname(basename='EMAT',iqmt=iq,filnam=fnemat)
     call genfilname(nodotpar=.true.,basename='X0_TIMING',bzsampl=bzsampl,&
          acont=acont,nar=.not.aresdf,iqmt=iq,procs=procs,rank=rank, &
          filnam=fnxtim)
     call genfilname(basename='X0',bzsampl=bzsampl,acont=acont,nar=.not.aresdf,&
          iqmt=iq,filnam=fnchi0)
     call genfilname(basename='X0',bzsampl=bzsampl,acont=acont,nar=.not.aresdf,&
          iqmt=iq,procs=procs,rank=rank,filnam=fnchi0_t)
  end if
  ! remove timing files from previous runs
  call filedel(trim(fnxtim))
  ! calculate k+q and G+k+q related variables
  call init1offs(qvkloff(1,iq))
  ! generate link array for tetrahedra
  if (tetradf) call gentetlinkp(vql(1,iq),tetraqweights)
  ! find highest (partially) occupied and lowest (partially) unoccupied states
  call findocclims(iq,istocc0,istocc,istunocc0,istunocc,isto0,isto,istu0,istu)
  ! find limits for band combinations
  call ematbdcmbs(emattype)
  ! check if q-point is Gamma point
  tq0=tqgamma(iq)
  if (tq0) then
     write(unitout,'(a)') 'Info('//trim(thisnam)//'): Gamma q-point: using &
          &momentum matrix elements for dielectric function'
  end if
  ! write out matrix size of response function
  write(unitout,'(a,i6)') 'Info('//thisnam//'): number of G+q vectors &
       &(local field effects):',ngq(iq)
  write(unitout,'(a,4i6)') 'Info('//thisnam//'): lowest (partially) &
       & unoccupied state: ',istunocc0
  write(unitout,'(a,4i6)') 'Info('//thisnam//'): highest (partially)&
       & occupied state  : ',istocc0
  write(unitout,'(a,4i6)') 'Info('//thisnam//'): limits for band combinations&
       & nst1,nst2,nst3,nst4:',nst1,nst2,nst3,nst4
  write(unitout,'(a,4i6)') 'Info('//thisnam//'): limits for band combinations&
       & istl1,istu1,istl2,istu2:',istl1,istu1,istl2,istu2
  write(unitout,'(a,4i6)') 'Info('//thisnam//'): limits for band combinations&
       & istl3,istu3,istl4,istu4:',istl3,istu3,istl4,istu4
  ! allocate arrays for eigenvalue and occupation number differences
  if (allocated(deou)) deallocate(deou)
  allocate(deou(nst1,nst2))
  if(allocated(deuo)) deallocate(deuo)
  allocate(deuo(nst3,nst4))
  if (allocated(docc12)) deallocate(docc12)
  allocate(docc12(nst1,nst2))
  if (allocated(docc21)) deallocate(docc21)
  allocate(docc21(nst3,nst4))
  ! allocate matrix elements arrays
  if (allocated(xiou)) deallocate(xiou)
  allocate(xiou(nst1,nst2,n))
  if (allocated(xiuo)) deallocate(xiuo)
  allocate(xiuo(nst3,nst4,n))
  if (allocated(pmou)) deallocate(pmou)
  allocate(pmou(3,nst1,nst2))
  if (allocated(pmuo)) deallocate(pmuo)
  allocate(pmuo(3,nst3,nst4))
  ! allocate arrays
  allocate(hdg(nst1,nst2,nkpt))
  allocate(scis12(nst1,nst2))
  allocate(scis21(nst2,nst1))
  allocate(w(nwdf))
  allocate(wreal(nwdfp))
  allocate(chi0h(3,3,nwdfp))
  allocate(chi0w(n,2,3,nwdfp))
  allocate(chi0(n,n,nwdfp))
  allocate(wou(nwdf))
  allocate(wuo(nwdf))
  allocate(wouw(nwdf),wuow(nwdf),wouh(nwdf),wuoh(nwdf))
  allocate(zvou(n),zvuo(n))
  allocate(bsedg(nst1,nst2))
  scis12(:,:)=0.d0
  scis21(:,:)=0.d0
  if (tetradf) then
     allocate(cw(nwdf),cwa(nwdf),cwsurf(nwdf))
     if (tetracw1k) allocate(cwt(nstsv,nstsv),cw1k(nst1,nst2,nwdfp), &
          cwa1k(nst1,nst2,nwdfp),cwsurf1k(nst1,nst2,nwdfp))
  end if
  ! generate complex energy grid
  call genwgrid(nwdf,wdos,acont,0.d0,w_cmplx=w)
  wreal(:)=dble(w(wi:wf))
  if (wreal(1).lt.epstetra) wreal(1)=epstetra
  ! initializations
  chi0(:,:,:)=zzero
  chi0w(:,:,:,:)=zzero
  chi0h(:,:,:)=zzero
  if (tscreen) then
     ! generate radial integrals wrt. sph. Bessel functions
     call ematrad(iq)
     ! delete timing information of previous runs
     call filedel(trim(fnetim))
     ! write information
     write(unitout,'(a,i6)') 'Info('//thisnam//'): number of G+q vectors:', &
          ngq(iq)
     call ematqalloc
  end if
  if (task.eq.345) then
     call getbsediag
     write(unitout,'("Info(",a,"): read diagonal of BSE kernel")') trim(thisnam)
     write(unitout,'(" mean value : ",2g18.10)') bsed
  end if
  bsedg(:,:)=bsed
  ! loop over k-points
  do ik=1,nkpt
     ! k-point analysis
     if (.not.transik(ik,dftrans)) cycle
     call chkpt(3,(/task,iq,ik/),'dfq: task, q-point index, k-point index')
     cpuosc=0.d0
     cpuupd=0.d0
     call cpu_time(cpu0)
     ikq=ikmapikq(ik,iq)
     call getdevaldoccsv(iq,ik,ikq,istl1,istu1,istl2,istu2,deou,docc12, &
          scis12)
     call getdevaldoccsv(iq,ik,ikq,istl2,istu2,istl1,istu1,deuo,docc21, &
          scis21)
     if (tscreen) then
        ! do not use scissors correction for screening
        if (task.eq.430) then
          scis12(:,:)=0.d0
          scis21(:,:)=0.d0
        end if
        ! for screening calculate matrix elements of plane wave on the fly
        call ematqk1(iq,ik)
        if (.not.allocated(xiuo)) allocate(xiuo(nst3,nst4,n))
        if (.not.allocated(pmuo)) allocate(pmuo(3,nst3,nst4))
     end if
     ! add BSE diagonal shift use with BSE-kernel
     if (task.eq.345) then
        scis12(:,:)=scis12(:,:)+bsedg(:,:)
        scis21(:,:)=scis21(:,:)+transpose(bsedg(:,:))
     end if
     ! get matrix elements (exp. expr. or momentum op.)
     call getpemat(iq,ik,trim(fnpmat),trim(fnemat),m12=xiou,m34=xiuo, &
          p12=pmou,p34=pmuo)
     ! set matrix elements to one for Lindhard function
     if (lindhard) then
       ! set G=0 components to one
       xiou(:,:,1)=zone
       xiuo(:,:,1)=zone
       ! set G/=0 components to zero
       if (n.gt.1) then
         xiou(:,:,2:)=zzero
         xiuo(:,:,2:)=zzero
       end if
       ! set momentum matrix elements to one
       pmou(:,:,:)=zone
       pmuo(:,:,:)=zone
     end if
     if (tetracw1k) then
        do iw=1,nwdfp
           call tetcwifc_1k(ik,nkpt,nstsv,evalsv,efermi,wreal(iw),2,cwt)
           cw1k(:,:,iw)=cwt(istl1:istu1,istl2:istu2)
           call tetcwifc_1k(ik,nkpt,nstsv,evalsv,efermi,-wreal(iw),2,cwt)
           cwa1k(:,:,iw)=cwt(istl1:istu1,istl2:istu2)
           call tetcwifc_1k(ik,nkpt,nstsv,evalsv,efermi,wreal(iw),4,cwt)
           cwsurf1k(:,:,iw)=cwt(istl1:istu1,istl2:istu2)
        end do
     end if
     if (tscreen) then
        ! we don't need anti-resonant parts here, assign them the same
        ! value as for resonant parts, resulting in a factor of two.
        do igq=1,n
           xiuo(:,:,igq)=transpose(xiou(:,:,igq))
        end do
        do j=1,3
           pmuo(j,:,:)=transpose(pmou(j,:,:))
        end do
        deuo(:,:)=transpose(deou(:,:))
        docc21(:,:)=transpose(docc12(:,:))
        scis21(:,:)=transpose(scis12(:,:))
     end if
     ! turn off antiresonant terms (type 2-1 band combiantions) for Kohn-Sham
     ! response function
     if ((.not.aresdf).and.(.not.tscreen)) then
        xiuo(:,:,:)=zzero
        pmuo(:,:,:)=zzero
     end if
     do ist1=1,istocc0-istunocc0+1
        do ist2=1,istocc0-istunocc0+1
           j=ist1+istunocc0-1
           ! set lower triangle of first block to zero
           if (ist1.gt.ist2) then
              xiou(j,ist2,:)=zzero
              pmou(:,j,ist2)=zzero
           end if
           ! set diagonal to zero (project out intraband contributions)
           if ((.not.intraband).and.(ist1.eq.ist2)) then
              xiou(j,ist2,:)=zzero
              pmou(:,j,ist2)=zzero
           end if
           ! set upper triangle of second block to zero
           ! also set diagonal to zero to avoid double counting
           if (ist1.ge.ist2) then
              xiuo(ist2,j,:)=zzero
              pmuo(:,ist2,j)=zzero
           end if
        end do
     end do
     call cpu_time(cpu1)
     cpuread=cpu1-cpu0
     do ist1=1,nst1
        do ist2=1,nst2
           !---------------------!
           !     denominator     !
           !---------------------!
           ! absolute band indices
           i1=ist1
           i2=istunocc0+ist2-1
           ! band analysis
           if (.not.transijst(ik,i1,i2,dftrans)) cycle
           call cpu_time(cpu0)
           ! user request termination
           call terminateqry('dfq')
           if (tetradf) then
              ! mirror index pair on diagonal if necessary
              if (i1.gt.i2) then
                 j1=ist2
                 j2=ist1-istunocc0+1
              else
                 j1=ist1
                 j2=ist2
              end if
              ! read weights for tetrahedron method
              if (tetracw1k) then
                 cw(wi:wf)=cw1k(ist1,ist2,:)
                 cwa(wi:wf)=cwa1k(ist1,ist2,:)
                 cwsurf(wi:wf)=cwsurf1k(ist1,ist2,:)
              else
                 call gettetcw(iq,ik,j1,j2,nst1,nst2,nwdf, &
                   trim(fnwtet),cw,cwa,cwsurf)
              end if
              ! include occupation number differences
              wou(wi:wf)=docc12(ist1,ist2)*cmplx(cw(wi:wf),cwsurf(wi:wf),8)/ &
                   omega
              wuo(wi:wf)=-docc21(ist2,ist1)*cmplx(cwa(wi:wf),0.d0,8)/omega
              if (tq0) then
                 ! rescale: use delta-function delta(e_nmk + scis_nmk - w)
                 wouw(wi:wf)=cmplx(dble(wou(wi:wf)),aimag(wou(wi:wf))* &
                      deou(ist1,ist2)/ &
                      (-wreal(:)-scis12(ist1,ist2)))
                 wuow(wi:wf)=cmplx(dble(wuo(iw:wf)),aimag(wuo(wi:wf))* &
                      deuo(ist2,ist1)/ &
                      (-wreal(:)-scis21(ist2,ist1)))
                 wouh(wi:wf)=cmplx(dble(wou(wi:wf)),aimag(wou(wi:wf))* &
                      deou(ist1,ist2)**2/ &
                      (-wreal(:)-scis12(ist1,ist2))**2)
                 wuoh(wi:wf)=cmplx(dble(wuo(wi:wf)),aimag(wuo(wi:wf))* &
                      deuo(ist2,ist1)**2/ &
                      (-wreal(:)-scis21(ist2,ist1))**2)
              end if
           else
              ! include occupation number differences
              wou(wi:wf)=docc12(ist1,ist2)*wkpt(ik)/omega/(w(wi:wf)+ &
                   deou(ist1,ist2)+scis12(ist1,ist2)+zi*brd)
              wuo(wi:wf)=docc21(ist2,ist1)*wkpt(ik)/omega/(w(wi:wf)+ &
                   deuo(ist2,ist1)+scis21(ist2,ist1)+tordf*zi*brd)
              wouw(wi:wf)=wou(wi:wf)
              wuow(wi:wf)=wuo(wi:wf)
              wouh(wi:wf)=wou(wi:wf)
              wuoh(wi:wf)=wuo(wi:wf)
           end if
           !----------------------------------!
           !     update response function     !
           !----------------------------------!
           zvou(:)=xiou(ist1,ist2,:)
           zvuo(:)=xiuo(ist2,ist1,:)
           do iw=wi,wf
              ! body
              call zgerc(n,n,wou(iw),zvou,1,zvou,1,chi0(:,:,iw-wi+1),n)
              call zgerc(n,n,wuo(iw),zvuo,1,zvuo,1,chi0(:,:,iw-wi+1),n)
              if (tq0) then
                 do oct1=1,3
                    ! wings
                    chi0w(2:,1,oct1,iw-wi+1)=chi0w(2:,1,oct1,iw-wi+1)+ &
                          wouw(iw)*pmou(oct1,ist1,ist2)*conjg(zvou(2:))+ &
                          wuow(iw)*pmuo(oct1,ist2,ist1)*conjg(zvuo(2:))
                    chi0w(2:,2,oct1,iw-wi+1)=chi0w(2:,2,oct1,iw-wi+1)+ &
                          wouw(iw)*zvou(2:)*conjg(pmou(oct1,ist1,ist2))+ &
                          wuow(iw)*zvuo(2:)*conjg(pmuo(oct1,ist2,ist1))
                    do oct2=1,3
                       ! head
                       chi0h(oct1,oct2,iw-wi+1)=chi0h(oct1,oct2,iw-wi+1)+ &
                            wouh(iw)*pmou(oct1,ist1,ist2)* &
                            conjg(pmou(oct2,ist1,ist2))+ &
                            wuoh(iw)*pmuo(oct1,ist2,ist1)* &
                            conjg(pmuo(oct2,ist2,ist1))
                    end do
                 end do
              end if
           end do
           call cpu_time(cpu0)
           cpuupd=cpuupd+cpu0-cpu1
           ! end loop over states combinations
        end do
     end do
     cputot=cpuread+cpuosc+cpuupd
     ! timing information
     call dftim(iq,ik,trim(fnxtim),cpuread,cpuosc,cpuupd, &
          cputot)
     ! synchronize
     if (.not.tscreen) call barrier
     ! end loop over k-points
  end do
  if (tscreen) call ematqdealloc
  ! symmetrize head
  if (tq0) then
     allocate(chi0hs(3,3,nwdfp))
     do oct1=1,3
        do oct2=1,3
           chi0hs(oct1,oct2,:)=zzero
           do i=1,3
              do j=1,3
                 chi0hs(oct1,oct2,:)=chi0hs(oct1,oct2,:)+symt2(oct1,oct2,i,j)* &
                      chi0h(i,j,:)
              end do
           end do
        end do
     end do
     chi0h(:,:,:)=chi0hs(:,:,:)
     deallocate(chi0hs)
  end if
  ! write dielectric tensor to file
  if (rank.eq.0) call writedielt(nwdf,dble(w),chi0h,0)
  ! write response function to file
  if (tscreen) then
     ! write out screening
     call getunit(un)
     open(un,file=trim(fnscreen),form='formatted',action='write', &
          status='replace')
     call putscreen(un,tq0,n,chi0h(:,:,1),chi0w(:,:,:,1),chi0(:,:,1))
     call writevars(un,iq,0)
     close(un)
  else
     do j=0,procs-1
        if (rank.eq.j) then
           do iw=wi,wf
              call putx0(tq0,iq,iw-wi+1,trim(fnchi0_t),'',&
                   chi0(:,:,iw-wi+1),chi0w(:,:,:,iw-wi+1),chi0h(:,:,iw-wi+1))
           end do
        end if
        call barrier
     end do
  end if

  deallocate(chi0,chi0h,chi0w)
  deallocate(docc12,docc21,scis12,scis21)
  deallocate(deou,deuo,wou,wuo,wouw,wuow,wouh,wuoh,zvou,zvuo)
  deallocate(xiou,xiuo,pmou,pmuo)
  deallocate(bsedg)
  deallocate(w,wreal)
  if (tetradf) then
     deallocate(cw,cwa,cwsurf)
     if (tetracw1k) deallocate(cwt,cw1k,cwa1k,cwsurf1k)
  end if
end subroutine dfq
!EOC
