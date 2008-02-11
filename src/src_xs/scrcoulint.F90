
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine scrcoulint
  use modmain
  use modmpi
  use modxs
  use m_genfilname
  use m_tdgauntgen
  use m_findgntn0
  use m_writegqpts
  implicit none
  ! local variables
  character(*), parameter :: thisnam='scrcoulint'
  real(8), parameter :: epsortho=1.d-12
  integer :: iknr,jknr,iqr,iq,isym,jsym,igq1,igq2,n,iflg,flg,j
  integer :: ngridkt(3),iv(3),ivgsym(3),ivg1(3),ivg2(3),lspl
  logical :: nosymt,reducekt,tq0
  real(8) :: vklofft(3),vqr(3),vq(3),vtl(3),v2(3),s(3,3),si(3,3),t1,t2,t3
  real(8), allocatable :: potcl(:,:,:)
  integer, allocatable :: igqmap(:,:),isyma(:,:),ivgsyma(:,:,:),nsyma(:)
  complex(8), allocatable :: phf(:,:,:)
  logical, allocatable :: done(:)
  real(8), external :: r3taxi
  logical, external :: tqgamma
  ! save global variables
  nosymt=nosym
  reducekt=reducek
  ngridkt(:)=ngridk(:)
  vklofft(:)=vkloff(:)
  ! map variables for screened Coulomb interaction
  call initbse
  nosym=nosymscr
  ! no symmetries implemented for screened Coulomb interaction
  reducek=.false.
  ! q-point set of screening corresponds to (k,kp)-pairs
  ngridk(:)=ngridq(:)
  vkloff(:)=vkloffbse(:)
  if (nemptyscr.eq.-1) nemptyscr=nempty
  emattype=2
  call init0
  call init1
  call init2xs
  ! read Fermi energy from file
  call readfermi
  ! save variables for the Gamma q-point
  call tdsave0
  ! generate Gaunt coefficients
  call tdgauntgen(lmaxapw,lmaxemat,lmaxapw)
  ! find indices for non-zero Gaunt coefficients
  call findgntn0(lmaxapwtd,lmaxapwtd,lmaxemat,tdgnt)
  write(unitout,'(a,3i8)') 'Info('//thisnam//'): Gaunt coefficients generated &
       &within lmax values:', lmaxapw,lmaxemat,lmaxapw
  write(unitout,'(a,i6)') 'Info('//thisnam//'): number of q-points: ',nqpt
  call flushifc(unitout)
  call genfilname(dotext='_SCR.OUT',setfilext=.true.)
  call findocclims(0,istocc0,istocc,istunocc0,istunocc,isto0,isto,istu0,istu)
  call ematbdcmbs(emattype)
  call genfilname(dotext='_SCI.OUT',setfilext=.true.)
  ! check number of empty states
  if (nemptyscr.lt.nempty) then
     write(*,*)
     write(*,'("Error(",a,"): too few empty states in screening eigenvector &
          &file - the screening should include many empty states &
          &(BSE/screening)",2i8)') trim(thisnam),nempty,nemptyscr
     write(*,*)
     call terminate
  end if
  if (rank.eq.0) then
     call writekpts
     call writeqpts
  end if

  flg=2
  allocate(done(nqpt))
  allocate(nsyma(nqpt),isyma(maxsymcrys,nqpt),ivgsyma(3,maxsymcrys,nqpt))
  allocate(igqmap(ngqmax,nqpt))
  done(:)=.false.
  ! loop over non-reduced number of k-points
  do iknr=1,nkptnr     
     do jknr=iknr,nkptnr
        iv(:)=ivknr(:,jknr)-ivknr(:,iknr)
        iv(:)=modulo(iv(:),ngridk(:))
        ! q-point (reduced)
        iqr=iqmapr(iv(1),iv(2),iv(3))
        vqr(:)=vqlr(:,iqr)
        ! q-point (non-reduced)
        iq=iqmap(iv(1),iv(2),iv(3))
        tq0=tqgamma(iq)
        vq(:)=vql(:,iq)

! first Brillouin zone q-point
call mapto1bz(vq,v2)
write(20,'(a,i6,3f12.3,3x,3f12.3,3x,2f12.3)') 'q,BZ,1BZ',iq,vq,v2, &
     sqrt(sum(matmul(bvec,vq)**2)),sqrt(sum(matmul(bvec,v2)**2))

        ! local field effects size
        n=ngq(iq)
        allocate(phf(nqpt,n,n),potcl(nqpt,n,n))

        ! symmetry that transforms non-reduced q-point to reduced one
        nsyma(iq)=0
        do isym=1,nsymcrys
           lspl=lsplsymc(isym)
           s(:,:)=dble(symlat(:,:,lspl))
           call r3mtv(s,vqr,v2)
           call r3frac(epslat,v2,ivgsym)
           t1=r3taxi(vq,v2)
           if (t1.lt.epslat) then
              nsyma(iq)=nsyma(iq)+1
              isyma(nsyma(iq),iq)=isym
              ivgsyma(:,nsyma(iq),iq)=-ivgsym(:)
           end if
        end do
        
!!$        do j=1,nsyma(iq)
!!$           isym=isyma(j,iq)
!!$           lspl=lsplsymc(isym)
!!$           do igq1=1,n
!!$              ivg1(:)=ivg(:,igqig(igq1,iq))
!!$              ! iv = sLT * ( G + G_s ); L...lattice; T...transpose
!!$              iv=matmul(transpose(symlat(:,:,lspl)),ivg1+ivgsyma(:,j,iq))
!!$              v2(:)=dble(iv(1))*bvec(:,1)+dble(iv(2))*bvec(:,2) &
!!$                   +dble(iv(3))*bvec(:,3)
!!$              t1=v2(1)**2+v2(2)**2+v2(3)**2
!!$write(*,'(a,4i6,g18.10)') 'reduce:',iq,j,isym,igq1,sqrt(t1)
!!$              if (t1.gt.gqmax**2) goto 10
!!$           end do
!!$           jsym=isym
!!$           ivgsym(:)=ivgsyma(:,j,iq)
!!$           s(:,:)=dble(symlat(:,:,lspl))
!!$           goto 20
!!$10         continue
!!$        end do
!!$        write(*,*)
!!$        write(*,'("Error(",a,"): failed to reduce q-point: ",i8)') &
!!$             trim(thisnam),iq
!!$        write(*,*)
!!$        call terminate
!!$20      continue


!!           igqmap(igq1,iq)=ivgigq(iv(1),iv(2),iv(3),iq)


write(*,'(a,3i6,3x,192i4)') 'ik1,ik2,iq,symops(iq)',iknr,jknr,iq, &
     isyma(1:nsyma(iq),iq)

        ! cross check symmetry relation (vq = a_isym vqr + G_isym)
        v2=matmul(transpose(symlat(:,:,lsplsymc(isyma(1,iq)))),vqr)+dble(ivgsyma(:,1,iq))
        v2=vq-v2
        if (any(abs(v2).gt.epslat)) then
write(*,'(a,2i6,3x,2i6,3f12.3)') 'deviation:',iknr,jknr,iqr,iq,v2
        end if

        call genfilname(iq=iq,dotext='_SCI.OUT',setfilext=.true.)
!!$        if (.not.done(iq)) call writegqpts(iq)
        call genfilname(dotext='_SCR.OUT',setfilext=.true.)
        ! calculate matrix elements
!!$        call init1xs(qvkloff(1,iq))
!!$        call ematrad(iq)
!!$        call ematqalloc
!!$        call ematqk1(iq,iknr)
!!$        call ematqdealloc
        call genfilname(dotext='_SCI.OUT',setfilext=.true.)

       

        if (.not.done(iq)) then
           ! calculate phase factor for dielectric matrix
           do igq1=1,n
              ivg1(:)=ivg(:,igqig(igq1,iq))
              ! find index transformation for G-vectors
              iv(:)=matmul(transpose(symlat(:,:,lsplsymc(scimap(isyma(1,iq))))),ivg1+ivgsym)
              igqmap(igq1,iq)=ivgigq(iv(1),iv(2),iv(3),iq)

write(*,'(a,5i6)') 'iknr,jknr,iq,igq1,igq1map',iknr,jknr,iq,igq1,igqmap(igq1,iq)

              do igq2=igq1,n
                 ! G-vector difference
                 ivg2(:)=ivg(:,igqig(igq2,iq))-ivg1(:)
                 ! translation vector s^-1*vtl(s^-1)
                 vtl=matmul(transpose(s),vtlsymc(:,scimap(isyma(1,iq))))
                 call r3frac(epslat,vtl,iv)
                 t1=twopi*dot_product(dble(ivg2),vtl)
                 t2=cos(t1)
                 t3=sin(t1)
                 if (abs(t2).lt.epsortho) t2=0.d0
                 if (abs(t3).lt.epsortho) t3=0.d0
                 ! phase factor
                 phf(iq,igq1,igq2)=cmplx(t2,t3,8)
                 phf(iq,igq2,igq1)=conjg(phf(iq,igq1,igq2))
!write(*,*) 'q,g,gp,phf',iq,igq1,igq2,phf(iq,igq1,igq2)
                 ! calculate weights for Coulomb potential
                 iflg=0
                 ! integrate weights for q=0 head and wings
                 ! and for q/=0 head
                 if (tq0) then
                    if (.not.((igq1.ne.1).and.(igq2.ne.1))) iflg=flg
                 end if
!!$                 call genwiq2xs(iflg,iq,igq1,igq2,potcl(iq,igq1,igq2))
                 potcl(iq,igq2,igq1)=potcl(iq,igq1,igq2)
if (iflg.ne.0) write(*,'(a,6i8,2g18.10)') 'ik,jk,q,flg,g,gp,potcl',iknr,jknr,iq,iflg,igq1,igq2,potcl(iq,igq1,igq2),fourpi/(gqc(igq1,iq)*gqc(igq2,iq))
              end do
           end do
           ! ** end if (.not.done(iq))
        end if



        




        ! *** read dielectric matrix and invert for >>iqr<<


        ! deallocate
        deallocate(phf,potcl)
        done(iq)=.true.
        ! end loops over non-reduced k-point combinations
     end do
  end do


write(*,*) 'minimum number of symmetry operation for q-point',minval(nsyma)
write(*,*) 'maximum number of symmetry operation for q-point',maxval(nsyma)



  call findgntn0_clear
  deallocate(done,isyma,nsyma,ivgsyma,igqmap)


  ! restore global variables
  nosym=nosymt
  reducek=reducekt
  ngridk(:)=ngridkt(:)
  vkloff(:)=vklofft(:)
  write(unitout,'(a)') "Info("//trim(thisnam)//"): Screening finished"
end subroutine scrcoulint



subroutine mapto1bz(vl,vl1bz)
  use modmain
  real(8), intent(in) :: vl(3)
  real(8), intent(out) :: vl1bz(3)
  ! local variables
  real(8) :: v0(3),v1(3)
  ! map the q-vector into the first Brillouin zone
  t1=1.d8
  do i1=-1,1
     do i2=-1,1
        do i3=-1,1
           v0=vl+dble((/i1,i2,i3/))
           v1=matmul(bvec,v0)
           t2=v1(1)**2+v1(2)**2+v1(3)**2
           ! favour positive values
           if (t2.lt.(t1+1.d-8)) then
              t1=t2
              vl1bz(:)=v0(:)
           end if
        end do
     end do
  end do
end subroutine mapto1bz
