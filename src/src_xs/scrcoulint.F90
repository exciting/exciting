
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine scrcoulint
  use modmain
  use modmpi
  use modxs
  use invert
  use m_tdgauntgen
  use m_findgntn0
  use m_writegqpts
  use m_genfilname
  use m_getunit

  implicit none
  ! local variables
  character(*), parameter :: thisnam='scrcoulint'
  real(8), parameter :: epsortho=1.d-12
  integer :: iknr,jknr,iqr,iq,iqrnr,isym,isymi,jsym,jsymi,igq1,igq2,n,iflg,flg,j
  integer :: ngridkt(3),iv(3),ivgsym(3),ivg1(3),ivg2(3),lspl,lspli,un
  integer :: idum1,idum2,idum3,oct,info
  logical :: nosymt,reducekt,tq0
  real(8) :: vklofft(3),vqr(3),vq(3),vtl(3),v2(3),s(3,3),si(3,3),t1,t2,t3
  real(8) :: rm(2,9)
  complex(8) :: scrnh0(3),scrnih0(3)
  character(256) :: fname
  real(8), allocatable :: potcl(:,:,:)
  complex(8), allocatable :: scrn(:,:),scrnw(:,:,:),scrnh(:)
  complex(8), allocatable :: scrni(:,:,:),tm(:,:),tmi(:,:)
  integer, allocatable :: igqmap(:,:),isyma(:,:),ivgsyma(:,:,:),nsyma(:)
  complex(8), allocatable :: phf(:,:,:)
  logical, allocatable :: done(:)
  real(8), external :: r3taxi
  integer, external :: octmap
  logical, external :: tqgamma

info=0

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

  ! *** read dielectric matrix and invert for >>iqr<<
  call genfilname(dotext='_SCR.OUT',setfilext=.true.)
  call getunit(un)
  allocate(scrni(ngqmax,ngqmax,nqptr))
  do iqr=1,nqptr
     ! locate reduced q-point in non-reduced set
     iv(:)=nint(vqlr(:,iqr)*ngridq(:))
     iqrnr=iqmap(iv(1),iv(2),iv(3))
     n=ngq(iqrnr)
     allocate(scrn(n,n),scrnw(n,2,3),scrnh(9),tm(n,n),tmi(n,n))
     tq0=tqgamma(iqrnr)
write(*,*) 'iqr,ngq',iqr,n
     call genfilname(basename='SCREEN',iq=iqr,filnam=fname)
     open(un,file=trim(fname),form='formatted',action='read',status='old')
     do igq1=1,n
        do igq2=1,n
           if (tq0) then
              if ((igq1.eq.1).and.(igq2.eq.1)) then
!                 read(un,*) (idum1,idum1,idum3,scrnh(j),j=1,9)
                 read(un,*) (idum1,idum1,idum3,rm(1,j),rm(2,j),j=1,9)
                 scrnh(:)=cmplx(rm(1,:),rm(2,:),8)
write(*,*) 'head*** read:',scrnh
              end if
              if ((igq1.eq.1).and.(igq2.ne.1)) then
!                 read(un,*) (idum1,idum2,idum3,scrnw(igq2,1,j),j=1,3)
                 read(un,*) (idum1,idum2,idum3,rm(1,j),rm(2,j),j=1,3)
                 scrnw(igq2,1,:)=cmplx(rm(1,:3),rm(2,:3),8)
              end if
              if ((igq1.ne.1).and.(igq2.eq.1)) then
!                 read(un,*) (idum1,idum2,idum3,scrnw(igq1,2,j),j=1,3)
                 read(un,*) (idum1,idum2,idum3,rm(1,j),rm(2,j),j=1,3)
                 scrnw(igq1,2,:)=cmplx(rm(1,:3),rm(2,:3),8)
              end if
              if ((igq1.ne.1).and.(igq2.ne.1)) read(un,*) idum1,idum2,idum3,&
!                   scrn(igq1,igq2)
                   rm(1,1),rm(2,1)
              scrn(igq1,igq2)=cmplx(rm(1,1),rm(2,1),8)
!!$              if ((igq1.eq.1).and.(igq2.eq.1)) read(un,*) &
!!$                      idum1,idum2,scrnh(:)
!!$              if ((igq1.eq.1).and.(igq2.ne.1)) read(un,*) &
!!$                   idum1,idum2,scrnw(igq2,1,:)
!!$              if ((igq1.ne.1).and.(igq2.eq.1)) read(un,*) &
!!$                   idum1,idum2,scrnw(igq1,2,:)
!!$              if ((igq1.ne.1).and.(igq2.ne.1)) read(un,*) &
!!$                   idum1,idum2,scrn(igq1,igq2)
           else
!              read(un,*) idum1,idum2,idum3,scrn(igq1,igq2)
              read(un,*) idum1,idum2,idum3,rm(1,1),rm(2,1)
              scrn(igq1,igq2)=cmplx(rm(1,1),rm(2,1),8)
           end if
        end do
     end do
     close(un)


    if ((iqr.eq.1).and.(n.ne.1)) then
       do igq1=1,n
          do igq2=1,n
             write(300+oct,'(2i8,2g18.10)') igq1,igq2,scrn(igq1,igq2)
          end do
       end do
    end if



     ! invert dielectric matrix for q-point going to zero in x-, y-, and
     ! z-direction for q=0
     if (tq0) then
        do oct=1,3
           ! index for diagonal tensor component
           j=octmap(oct,oct)
           scrn(1,1)=scrnh(j)
           ! keep head of dielectric matrix for q=0
           scrnh0(oct)=scrn(1,1)
           if (n.gt.1) then
              scrn(1,2:n)=scrnw(2:n,1,oct)
              scrn(2:n,1)=scrnw(2:n,2,oct)
           end if
           tm=scrn

write(*,*) 'Scr. Coul. Int.: head of diel. matrix, oct:',oct,scrnh0(oct)

!!$do igq1=1,n
!!$do igq2=igq1,n
!!$scrn(igq2,igq1)=scrn(igq1,igq2)
!!$tm(igq1,igq2)=scrn(igq1,igq2)
!!$tm(igq2,igq1)=scrn(igq1,igq2)
!!$end do
!!$end do

! set to unity matrix for input to zposv
tmi(:,:)=0.d0
forall(j=1:n) tmi(j,j)=1.d0




           

           ! call zposv('u',numgdm(i)+2,numgdm(i)+2,aa,maxgv+2,bb,maxgv+2,info)
           ! solve linear equation system for positive definite symmetric
           ! matrix

write(200,*) scrn

!!!           call zposv('u',n,n,tm,n,tmi,n,info)

call zinvert_lapack(tm,tmi)


           if (info.ne.0) then
              write(*,*)
              write(*,'("Error(",a,"): zposv returned non-zero info : ",I8)') &
                   trim(thisnam),info
              write(*,*)
              call terminate
           end if
           scrni(:,:,iqr)=tmi(:,:)
           ! keep head of inverse dielectric matrix for q=0
           scrnih0(oct)=scrni(1,1,iqr)
write(*,*) 'Scr. Coul. Int.: head of inv. diel. matrix, oct:',oct,scrnih0(oct)
        end do
        ! symmetrize inverse dielectric matrix wrt. directions in which q goes
        ! to zero
        call symsci0(1,scrnh0,scrnih0,scrni(1,1,iqr))
write(*,*) 'Scr. Coul. Int.: symm. head of inv. diel. matrix:',scrni(1,1,iqr)
     else

do igq1=1,n
do igq2=igq1,n
scrn(igq2,igq1)=scrn(igq1,igq2)
end do
end do

tm(:,:)=scrn(:,:)

! set to unity matrix for input to zposv
tmi(:,:)=0.d0
forall(j=1:n) tmi(j,j)=1.d0

        call zposv('u',n,n,tm,n,tmi,n,info)
        if (info.ne.0) then
           write(*,*)
           write(*,'("Error(",a,"): zposv returned non-zero info : ",I8)') &
                trim(thisnam),info
           write(*,*)
           call terminate
        end if
     end if
     scrni(:,:,iqr)=tmi(:,:)
     deallocate(scrn,scrnw,scrnh,tm,tmi)
     ! end loop over reduced q-points
  end do
  call genfilname(dotext='_SCI.OUT',setfilext=.true.)

  ! flag for integrating the singular terms in the screened Coulomb interaction
  flg=1
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
        ! locate reduced q-point in non-reduced set
        iv(:)=nint(vqr(:)*ngridq(:))
        iqrnr=iqmap(iv(1),iv(2),iv(3))

        ! local field effects size
        n=ngq(iq)
        allocate(phf(nqpt,n,n),potcl(n,n,nqpt))

        ! symmetries that transform non-reduced q-point to reduced one
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
        
        ! find map from G-vectors to rotated G-vectors
        do j=1,nsyma(iq)
           isym=isyma(j,iq)
           lspl=lsplsymc(isym)
           isymi=scimap(isym)
           lspli=lsplsymc(isymi)
           do igq1=1,n
              ivg1(:)=ivg(:,igqig(igq1,iq))
              ! G1 = s^-1 * ( G + G_s )
              iv=matmul(transpose(symlat(:,:,lspli)),ivg1+ivgsyma(:,j,iq))
              v2=matmul(bvec,iv+vqr)
              t1=sqrt(sum(v2**2))
!!$write(*,'(a,5i6,3g18.10,3x,3g18.10)') 'reduce:',iq,j,isym,&
!!$     igq1,ivgigq(iv(1),iv(2),iv(3),iqrnr),vgql(:,igq1,iq)- &
!!$     matmul(transpose(symlat(:,:,lspl)),vqr+iv)
              if (t1.gt.gqmax) then
                 write(*,*) '*** need one more symmetry operation'
                 goto 10
              end if
              ! map G1+q
              igqmap(igq1,iq)=ivgigq(iv(1),iv(2),iv(3),iqrnr)
              if (igqmap(igq1,iq).le.0) then
                 write(*,*)
                 write(*,'("Error(",a,"): failed to map rotated G-vector")') &
                      trim(thisnam)
                 write(*,'(" non-reduced q-point                    :",i8)') iq
                 write(*,'(" reduced q-point                        :",i8)') iqr
                 write(*,'(" reduced q-point in non-reduced set     :",i8)') &
                      iqrnr
                 write(*,'(" G+q-vector index (non-reduced q-point) :",i8)') &
                      igq1
                 write(*,'(" rotated G-vector                       :",3i8)') iv
                 write(*,*)
                 call terminate
              end if
              ! end loop over G+
           end do
           ivgsym(:)=ivgsyma(:,j,iq)
           jsym=isym
           s(:,:)=dble(symlat(:,:,lspl))
           si(:,:)=dble(symlat(:,:,lspli))
           goto 20
10         continue
        end do
        write(*,*)
        write(*,'("Error(",a,"): failed to reduce q-point: ",i8)') &
             trim(thisnam),iq
        write(*,*)
        call terminate
20      continue

        ! cross check symmetry relation (q1 = s^-1 * vqr + G_s)
        if (sum(abs(vq-(matmul(transpose(s),vqr)+dble(ivgsym)))).gt.epslat) then
           write(*,*) 'deviation:',iknr,jknr,iqr,iq,v2
           write(*,*)
           write(*,'("Error(",a,"): cross checking of symmetry reduction for &
                &q-vector failed:")') trim(thisnam)
           write(*,'(" non-reduced q-point                    :",i8)') iq
           write(*,'(" reduced q-point                        :",i8)') iqr
           write(*,'(" reduced q-point in non-reduced set     :",i8)') &
                iqrnr
           write(*,'(" crystal symmetry operation             :",i8)') isym
           write(*,'(" umklapp G-vector                       :",3i8)') ivgsym
           write(*,*)
           call terminate
        end if

        !**********************************************************************
        write(123,'(a,5i6,3x,i6,3x,3i5)') 'iknr,jknr,iq,iqr,iqrnr,isym,ivgsym',&
             iknr,jknr,iq,iqr,iqrnr,isym,ivgsym
        !**********************************************************************

        if (.not.done(iq)) then
           do igq1=1,n
              ivg1(:)=ivg(:,igqig(igq1,iq))
!!$write(*,'(a,5i6)') 'iknr,jknr,iq,igq1,igq1map',iknr,jknr,iq,igq1,&
!!$     igqmap(igq1,iq)
              do igq2=igq1,n
                 ! G-vector difference
                 ivg2(:)=ivg(:,igqig(igq2,iq))-ivg1(:)
                 ! translation vector s^-1*vtl(s^-1)
                 vtl=matmul(transpose(s),vtlsymc(:,isymi))
                 call r3frac(epslat,vtl,iv)
                 t1=twopi*dot_product(dble(ivg2),vtl)
                 t2=cos(t1)
                 t3=sin(t1)
                 if (abs(t2).lt.epsortho) t2=0.d0
                 if (abs(t3).lt.epsortho) t3=0.d0
                 ! phase factor for dielectric matrix
                 phf(iq,igq1,igq2)=cmplx(t2,t3,8)
                 phf(iq,igq2,igq1)=conjg(phf(iq,igq1,igq2))
write(40,'(a,i5,2x,2i5,2x,2i5,2g18.10)') 'q,g,gp,isym,isymi,phf',iq,igq1,igq2,isym,isymi, &
     phf(iq,igq1,igq2)
                 ! calculate weights for Coulomb potential
                 iflg=0
                 ! integrate weights for q=0 for the head and wings
                 ! and for q/=0 for the head
                 if (tq0) then
                    if (.not.((igq1.ne.1).and.(igq2.ne.1))) iflg=flg
                 end if
                 call genwiq2xs(iflg,iq,igq1,igq2,potcl(igq1,igq2,iq))
                 potcl(igq2,igq1,iq)=potcl(igq1,igq2,iq)
if (iflg.ne.0) &
     write(*,'(a,6i8,2g18.10)') 'ik,jk,q,flg,g,gp,potcl',iknr,jknr,iq,iflg,igq1,igq2,potcl(igq1,igq2,iq),fourpi/(gqc(igq1,iq)*gqc(igq2,iq))
              end do
           end do
        end if

        call genfilname(iq=iq,dotext='_SCI.OUT',setfilext=.true.)
        if (.not.done(iq)) call writegqpts(iq)
        call genfilname(dotext='_SCR.OUT',setfilext=.true.)

        ! calculate matrix elements
!!$        call init1xs(qvkloff(1,iq))
!!$        call ematrad(iq)
!!$        call ematqalloc
!!$        call ematqk1(iq,iknr)
!!$        call ematqdealloc
        call genfilname(dotext='_SCI.OUT',setfilext=.true.)





        ! deallocate
        deallocate(phf,potcl)
        done(iq)=.true.
        ! end loops over non-reduced k-point combinations
     end do
  end do


  write(*,*) 'minimum number of symmetry operation for q-point',minval(nsyma)
  write(*,*) 'maximum number of symmetry operation for q-point',maxval(nsyma)

  call findgntn0_clear
  deallocate(done,isyma,nsyma,ivgsyma,igqmap,scrni)

  ! restore global variables
  nosym=nosymt
  reducek=reducekt
  ngridk(:)=ngridkt(:)
  vkloff(:)=vklofft(:)
  write(unitout,'(a)') "Info("//trim(thisnam)//"): Screening finished"
end subroutine scrcoulint



subroutine mapto1bz(vl,vl1bz,iv)
  use modmain
  real(8), intent(in) :: vl(3)
  real(8), intent(out) :: vl1bz(3)
  integer, intent(out) :: iv(3)
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
           ! favour positive coordinates
           if (t2.lt.(t1+1.d-8)) then
              t1=t2
              vl1bz(:)=v0(:)
              iv=-(/i1,i2,i3/)
           end if
        end do
     end do
  end do
end subroutine mapto1bz
