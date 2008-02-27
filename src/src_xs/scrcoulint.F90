
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
  integer :: idum1,idum2,idum3,oct1,oct,info
  logical :: nosymt,reducekt,tq0
  real(8) :: vklofft(3),vqr(3),vq(3),vtl(3),v2(3),s(3,3),si(3,3),t1,t2,t3
  real(8) :: rm(2,9)
  complex(8) :: scrnh0(3),scrnih0(3)
  character(256) :: fname
  real(8), allocatable :: potcl(:,:,:)
  complex(8), allocatable :: scrn(:,:),scrnw(:,:,:),scrnh(:)
  complex(8), allocatable :: scrni(:,:,:),tm(:,:),tmi(:,:)
  complex(8), allocatable :: phf(:,:,:)
  integer, allocatable :: igqmap(:,:),isyma(:,:),ivgsyma(:,:,:),nsyma(:)
  logical, allocatable :: done(:)
  real(8), external :: r3taxi
  integer, external :: octmap
  logical, external :: tqgamma

integer :: lmax1,lmax2,lmax3
real(8) :: cpu0,cpu1,cpu2,cpu3
real(8) :: cpu_init1xs,cpu_ematrad,cpu_ematqalloc,cpu_ematqk1,cpu_ematqdealloc

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

  ! read dielectric matrix and invert for reduced q-point set
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
     ! read in screening
     call genfilname(basename='SCREEN',iq=iqr,filnam=fname)
     open(un,file=trim(fname),form='formatted',action='read',status='old')
     do igq1=1,n
        do igq2=1,n
           if (tq0) then
              if ((igq1.eq.1).and.(igq2.eq.1)) then
                 read(un,*) (idum1,idum2,idum3,rm(1,j),rm(2,j),j=1,9)
                 scrnh(:)=cmplx(rm(1,:),rm(2,:),8)
              end if
              if ((igq1.eq.1).and.(igq2.ne.1)) then
                 read(un,*) (idum1,idum2,idum3,rm(1,j),rm(2,j),j=1,3)
                 scrnw(igq2,1,:)=cmplx(rm(1,:3),rm(2,:3),8)
              end if
              if ((igq1.ne.1).and.(igq2.eq.1)) then
                 read(un,*) (idum1,idum2,idum3,rm(1,j),rm(2,j),j=1,3)
                 scrnw(igq1,2,:)=cmplx(rm(1,:3),rm(2,:3),8)
              end if
              if ((igq1.ne.1).and.(igq2.ne.1)) read(un,*) idum1,idum2,idum3,&
                   rm(1,1),rm(2,1)
              scrn(igq1,igq2)=cmplx(rm(1,1),rm(2,1),8)
           else
              read(un,*) idum1,idum2,idum3,rm(1,1),rm(2,1)
              scrn(igq1,igq2)=cmplx(rm(1,1),rm(2,1),8)
           end if
        end do
     end do
     close(un)
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
           ! store screening in temporary array
           tm=scrn
           write(*,'(a,i5,2g18.10)') 'optcomp,    1/eps_00(q=0)      :', &
                oct,1.d0/scrnh0(oct)
           write(*,'(a,i5,2g18.10)') 'optcomp,      eps_00(q=0)      :', &
                oct,scrnh0(oct)
           call zinvert_hermitian(scrherm,tm,tmi)
           scrni(:,:,iqr)=tmi(:,:)
           ! keep head of inverse dielectric matrix for q=0
           scrnih0(oct)=scrni(1,1,iqr)
           write(*,'(a,i5,2g18.10)') 'optcomp,   eps^-1_00(q=0)      :', &
                oct,scrnih0(oct)
           write(*,'(a,i5,2g18.10)') 'optcomp, 1/eps^-1_00(q=0)      :', &
                oct,1.d0/scrnih0(oct)
        end do
        ! symmetrize head of inverse dielectric matrix wrt. the directions
        ! in which q goes to zero
        call symsci0(bsediagsym,scrnh0,scrnih0,scrni(1,1,iqr))
        write(*,'(a,i5,2g18.10)') 'optcomp, symm.   eps^-1_00(q=0):', &
             iqr,scrni(1,1,iqr)
        write(*,'(a,i5,2g18.10)') 'optcomp, symm. 1/eps^-1_00(q=0):', &
             iqr,1.d0/scrni(1,1,iqr)
     else
        tm(:,:)=scrn(:,:)
        write(*,'(a,i5,2g18.10)') 'iq,    1/eps_00(q)             :', &
             iqr,1.d0/scrn(1,1)
        write(*,'(a,i5,2g18.10)') 'iq,      eps_00(q)             :', &
             iqr,scrn(1,1)
        call zinvert_hermitian(scrherm,tm,tmi)
        scrni(:,:,iqr)=tmi(:,:)
     end if
     write(*,'(a,i5,2g18.10)') 'diel. matr.: iq,   eps^-1_00(q):', &
          iqr,scrni(1,1,iqr)
     write(*,'(a,i5,2g18.10)') 'diel. matr.: iq, 1/eps^-1_00(q):', &
          iqr,1.d0/scrni(1,1,iqr)
     write(*,*)
     deallocate(scrn,scrnw,scrnh,tm,tmi)
     ! end loop over reduced q-points
  end do

!!$  do iq=1,nqpt
!!$write(*,*) 'radial integrals for q-point:',iq
!!$     ! calculate radial integrals
!!$     call ematrad(iq )
!!$     call genfilname(basename='EMATRAD',iq=iq,filnam=fname)
!!$     call getunit(un)
!!$     open(un,file=trim(fname),form='unformatted',action='write', &
!!$          status='replace')
!!$     write(un) riaa,riloa,rilolo
!!$     close(un)
!!$  end do

  call genfilname(dotext='_SCI.OUT',setfilext=.true.)
  ! flag for integrating the singular terms in the screened Coulomb interaction
  flg=bsediagweight
  allocate(done(nqpt))
  allocate(nsyma(nqpt),isyma(maxsymcrys,nqpt),ivgsyma(3,maxsymcrys,nqpt))
  allocate(igqmap(ngqmax,nqpt))
  done(:)=.false.
  ! loop over non-reduced number of k-points
  do iknr=1,nkptnr
     do jknr=iknr,nkptnr


call cpu_time(cpu2)
cpu_init1xs=0.d0
cpu_ematrad=0.d0
cpu_ematqalloc=0.d0
cpu_ematqk1=0.d0
cpu_ematqdealloc=0.d0



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

        ! symmetries that transform non-reduced q-point to reduced one, namely
        ! q1 = s^-1 * q + G_s. Here, q1 is vq, q is vqr.
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
              ! |G1 + q|
              v2=matmul(bvec,iv+vqr)
              t1=sqrt(sum(v2**2))
!!$write(*,'(a,5i6,3g18.10,3x,3g18.10)') 'reduce:',iq,j,isym,&
!!$     igq1,ivgigq(iv(1),iv(2),iv(3),iqrnr),vgql(:,igq1,iq)- &
!!$     matmul(transpose(symlat(:,:,lspl)),vqr+iv)
              if (t1.gt.gqmax) then
                 write(*,*) '*** need one more symmetry operation'
                 goto 10
              end if
              ! locate G1 + q in G+q-vector set
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
           ! store G1 vector
           ivgsym(:)=ivgsyma(:,j,iq)
           jsym=isym
           ! store symmetry
           s(:,:)=dble(symlat(:,:,lspl))
           ! store inverse of symmetry
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
        ! cross check symmetry relation (q1 = s^-1 * q + G_s)
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

        ! set up Coulomb potential and phase factor
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
                 ! phase factor for dielectric matrix (due to translations)
                 phf(iq,igq1,igq2)=cmplx(t2,t3,8)
                 phf(iq,igq2,igq1)=conjg(phf(iq,igq1,igq2))
write(40,'(a,i5,2x,2i5,2x,2i5,2g18.10)') 'q,g,gp,isym,isymi,phf',iq,igq1,igq2,isym,isymi,phf(iq,igq1,igq2)
                 ! calculate weights for Coulomb potential
                 iflg=0
                 ! integrate weights for q=0 for the head and wings
                 ! (and for q/=0 for the head ?? good idea??)
                 if (tq0) then
                    if (.not.((igq1.ne.1).and.(igq2.ne.1))) iflg=flg
                 end if
                 call genwiq2xs(iflg,iq,igq1,igq2,potcl(igq1,igq2,iq))
                 potcl(igq2,igq1,iq)=potcl(igq1,igq2,iq)
if (iflg.ne.0) &
write(*,'(a,6i8,2g18.10)') 'ik,jk,q,flg,g,gp,potcl',iknr,jknr,iq,iflg,igq1,igq2,potcl(igq1,igq2,iq),fourpi/(gqc(igq1,iq)*gqc(igq2,iq))
                 ! end loop over (G,Gp)-vectors
              end do
           end do
        end if

        call genfilname(iq=iq,dotext='_SCI.OUT',setfilext=.true.)
        if (.not.done(iq)) call writegqpts(iq)
        call genfilname(dotext='_SCR.OUT',setfilext=.true.)

        ! calculate matrix elements of the plane wave
        call cpu_time(cpu0)
        call init1xs(qvkloff(1,iq))
        call cpu_time(cpu1)
        cpu_init1xs=cpu_init1xs+cpu1-cpu0

write(*,*) 'iknr,jknr,iq,ngq(iq)',iknr,jknr,iq,ngq(iq)

        call cpu_time(cpu0)
        ! get radial integrals
        call ematrad(iq)
!!!        call getematrad(iq)


        call cpu_time(cpu1)
        cpu_ematrad=cpu_ematrad+cpu1-cpu0

        call cpu_time(cpu0)
        call ematqalloc
        call cpu_time(cpu1)
        cpu_ematqalloc=cpu_ematqalloc+cpu1-cpu0

        call cpu_time(cpu0)
        call ematqk1(iq,iknr)
        call cpu_time(cpu1)
        cpu_ematqk1=cpu_ematqk1+cpu1-cpu0

        call cpu_time(cpu0)
        call ematqdealloc
        call cpu_time(cpu1)
        cpu_ematqdealloc=cpu_ematqdealloc+cpu1-cpu0

        call genfilname(dotext='_SCI.OUT',setfilext=.true.)

        ! deallocate
        deallocate(phf,potcl)
        done(iq)=.true.
        ! end loops over non-reduced k-point combinations


! check matrix elements
if ((iknr.eq.2).and.(jknr.eq.7)) then

write(6000,*) xiou
write(6001,*) xiuo
stop 'stopped here'

end if



call cpu_time(cpu3)
t3=cpu_ematqdealloc+cpu_ematqk1+cpu_ematqk1+cpu_ematrad+cpu_init1xs
write(*,'(a,f12.3)') 'init1xs     :',cpu_init1xs
write(*,'(a,f12.3)') 'ematrad     :',cpu_ematrad
write(*,'(a,f12.3)') 'ematqalloc  :',cpu_ematqalloc
write(*,'(a,f12.3)') 'ematqk1     :',cpu_ematqk1
write(*,'(a,f12.3)') 'ematqdealloc:',cpu_ematqdealloc
write(*,'(a,f12.3)') '*** sum     :',t3
write(*,'(a,f12.3)') '*** rest    :',cpu3-cpu2-t3
write(*,'(a,f12.3)') '*** overall :',cpu3-cpu2
write(*,*)

     ! end loop over (k,kp) pairs
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


subroutine getematrad(iq)
  use modmain
  use modxs
  use m_genfilname
  use m_getunit
  implicit none
  ! arguments
  integer, intent(in) :: iq
  ! local variables
  integer :: lmax1,lmax2,lmax3,un
  character(256) :: fname
  lmax1=lmaxapwtd
  lmax2=lmaxemat
  ! lmax1 and lmax3 should be the same!
  lmax3=lmax1
  if (allocated(riaa)) deallocate(riaa)
  allocate(riaa(0:lmax1,apwordmax,0:lmax3,apwordmax,0:lmax2,natmtot, &
       ngq(iq)))
  if (allocated(riloa)) deallocate(riloa)
  allocate(riloa(nlomax,0:lmax3,apwordmax,0:lmax2,natmtot,ngq(iq)))
  if (allocated(rilolo)) deallocate(rilolo)
  allocate(rilolo(nlomax,nlomax,0:lmax2,natmtot,ngq(iq)))
  call genfilname(basename='EMATRAD',iq=iq,filnam=fname)
  call getunit(un)
  open(un,file=trim(fname),form='unformatted',action='read', &
       status='old')
  read(un) riaa,riloa,rilolo
  close(un)
end subroutine getematrad
