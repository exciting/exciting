
! Copyright (C) 2002-2005 S. Sharma, J. K. Dewhurst and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine linopt
  use modmain
#ifdef TETRA
  use modtetra
#endif
#ifdef XS
  use modxs, only: emattype,xiou,qvkloff,istocc0,istocc,istunocc0,istunocc
  use modxs, only: isto0,isto,istu0,istu,evalsv0,occsv0,ngq
  use m_getemat
  use m_genfilname
#endif
  implicit none
  ! local variables
  integer ik,nsk(3),iw,jw
  integer isym,lspl,iop
  integer i,i1,i2,n,recl
  real(8), parameter :: eps=1.d-8
  real(8) wd(2),wplas,t1,t2
  character(256) fname
  ! allocatable arrays
  real(8), allocatable :: w(:)
  real(8), allocatable :: fw(:)
  real(8), allocatable :: g(:)
  real(8), allocatable :: cf(:,:)
  real(8), allocatable :: e(:,:)
  real(8), allocatable :: f(:,:)
  real(8), allocatable :: sc(:,:,:)
  real(8), allocatable :: d(:)
  real(8), allocatable :: eps1(:)
  real(8), allocatable :: eps2(:)
  real(8), allocatable :: sigma1(:)
  real(8), allocatable :: sigma2(:)
  real(8), allocatable :: delta(:,:)
  real(8), allocatable :: pmatint(:,:)
  complex(8), allocatable :: evecfv(:,:)
  complex(8), allocatable :: evecsv(:,:)
  complex(8), allocatable :: apwalm(:,:,:,:)
#ifdef XS
  integer, parameter :: iq=1
  logical :: tqfmt
#endif
  complex(8), allocatable :: pmat(:,:,:)
#ifdef TETRA
  integer :: m,ist1,ist2,iqnr,iv(3)
  real(8), parameter :: epstetra=1.d-8
  real(8) :: escal,sum1,sum2,vr(3)
  real(8), allocatable :: e1(:,:)
  real(8), allocatable :: cwsurf(:,:,:),cw(:,:,:),cwa(:,:,:)
#endif
  !<sag>
  real(8) :: t3,scis
  real(8), allocatable :: eps1r(:)
  logical :: tev
  integer :: bzsmpl
  logical, external :: tqgamma
  ! default sampling of Brillouin zone (trilinear method)
  bzsmpl=0
  ! output in electron volt
  tev=.true.
  optltz=(optswidth.ne.0.d0)
  if (optltz) bzsmpl=2
  !</sag>
#ifdef TETRA
  if (tetra.and.optltz) then
     write(*,*)
     write(*,'("Error(linopt): specified tetrahedron method and Lorentzian &
          &broadening")')
     write(*,*)
     stop
  end if
  if (tetra) bzsmpl=1
#endif
  if ((usegdft).and.(xctype.lt.0)) then
     write(*,*)
     write(*,'("Error(linopt): generalised DFT cannot be used with exact &
          &exchange")')
     write(*,*)
     stop
  end if
  ! initialise universal variables
  call init0
  call init1
emattype=0
  call init2xs
#ifdef XS
  tqfmt=.not.tqgamma(iq)
  if (tqfmt) then
     call tdsave0
     call genfilname(iq=iq,setfilext=.true.)
     ! take first q-point
     call init1xs(qvkloff(1,iq))
     if (allocated(evalsv0)) deallocate(evalsv0)
     allocate(evalsv0(nstsv,nkpt))
     if (allocated(occsv0)) deallocate(occsv0)
     allocate(occsv0(nstsv,nkpt))
     if (allocated(xiou)) deallocate(xiou)
     allocate(xiou(nstsv,nstsv,ngq(iq)))
     call findocclims(iq,istocc0,istocc,istunocc0,istunocc,isto0,isto,istu0, &
          istu)
     call ematbdcmbs(emattype)
  end if
#endif
  ! read Fermi energy from file
  call readfermi
  ! allocate local arrays
  allocate(w(nwdos))
  allocate(fw(nwdos))
  allocate(g(nwdos))
  allocate(cf(3,nwdos))
  n=nstsv*nstsv
  allocate(e(n,nkpt))
  allocate(f(n,nkpt))
  allocate(sc(3,3,nsymcrys))
  allocate(d(nsymcrys))
  allocate(eps1(nwdos),eps2(nwdos))
  allocate(sigma1(nwdos),sigma2(nwdos))
  ! allocate first-variational eigenvector array
  allocate(evecfv(nmatmax,nstfv))
  ! allocate second-variational eigenvector array
  allocate(evecsv(nstsv,nstsv))
  ! allocate the momentum matrix elements array
  allocate(pmat(3,nstsv,nstsv))
  ! allocate intraband matrix elements
  allocate(pmatint(nstsv,nkpt))
  ! set up for generalised DFT correction if required
  allocate(delta(nstsv,nstsv))
  if (bzsmpl.ge.1) then
     allocate(eps1r(nwdos))
     eps1r(:)=0.d0
  end if
#ifdef TETRA
  if (tetra) then
     allocate(e1(nstsv,nkpt))
     allocate(cw(nstsv,nstsv,nkpt))
     allocate(cwa(nstsv,nstsv,nkpt))
     allocate(cwsurf(nstsv,nstsv,nkpt))
  end if
#endif
  if (usegdft) then
     ! initialisation required for generalised DFT
     allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
     call readstate
     call poteff
     call linengy
     call genapwfr
     call genlofr
  end if
  ! pre-calculation for symmetrisation
  do isym=1,nsymcrys
     lspl=lsplsymc(isym)
     ! symmetry matrix in Cartesian coordinates
     sc(:,:,isym)=symlatc(:,:,lspl)
     !<sag> calculate 3,3 minor of symmetry element </sag>
     d(isym)=sc(1,1,isym)*sc(2,2,isym)-sc(1,2,isym)*sc(2,1,isym)
  end do
  ! energy interval should start from zero
  wdos(1)=0.d0
  ! generate energy grid
  t1=(wdos(2)-wdos(1))/dble(nwdos)
  do iw=1,nwdos
     w(iw)=t1*dble(iw-1)+wdos(1)
  end do
  ! number of subdivisions used for interpolation
  do i=1,3
     nsk(i)=max(ngrdos/ngridk(i),1)
  end do
  ! find the record length for momentum matrix element file
  inquire(iolength=recl) pmat
#ifdef XS
  if (.not.tqfmt)  then
#endif
     open(50,file='PMAT.OUT',action='READ',form='UNFORMATTED',access='DIRECT', &
          recl=recl)
#ifdef XS
  end if
#endif
  ! loop over number of desired optical components
  do iop=1,noptcomp
     i1=optcomp(1,iop)
     i2=optcomp(2,iop)
     ! open files for writting
     select case(bzsmpl)
     case(0)
        write(fname,'("EPSILON_",2I1,".OUT")') i1,i2
     case(1)
        write(fname,'("EPSILON_TET_",2I1,".OUT")') i1,i2
     case(2)
        write(fname,'("EPSILON_LTZ_",2I1,".OUT")') i1,i2
     end select
     open(60,file=trim(fname),action='WRITE',form='FORMATTED')
     select case(bzsmpl)
     case(0)
        write(fname,'("SIGMA_",2I1,".OUT")') i1,i2
     case(1)
        write(fname,'("SIGMA_TET_",2I1,".OUT")') i1,i2
     case(2)
        write(fname,'("SIGMA_LTZ_",2I1,".OUT")') i1,i2
     end select
     open(61,file=trim(fname),action='WRITE',form='FORMATTED')
     e(:,:)=0.d0
     f(:,:)=0.d0
     do ik=1,nkpt
        ! compute generalised DFT correction if required
        if (usegdft) then
           call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
           call getevecsv(vkl(1,ik),evecsv)
           call match(ngk(ik,1),gkc(1,ik,1),tpgkc(1,1,ik,1),sfacgk(1,1,ik,1),&
                apwalm)
           call gdft(ik,apwalm,evecfv,evecsv,delta)
        end if
        ! read matrix elements from direct-access file
#ifdef XS
        if (tqfmt) then
           call getemat(1,ik,.true.,'EMAT_Q00001.OUT',xiou)
           call linoptkpq(iq,ik,xiou,e(1,ik),f(1,ik))
        else
#endif
           read(50,rec=ik) pmat
           call linoptk(ik,i1,i2,sc,d,delta,pmat,e(1,ik),f(1,ik),pmatint(1,ik))
#ifdef XS
        end if
#endif
     end do
#if TETRA
     if (tetra) then
        if (tqfmt) then
!!$           ! get index to reducible q-point which is commensurate to k-point set
!!$           vr(:)=vql(:,iq)*ngridk(:)
!!$           call r3frac(epslat,vr,iv)
!!$           if (sum(abs(vr)).gt.epscomm) then
!!$              write(*,*)
!!$              write(*,'("Error(linopt): q-point not commensurate with k-point &
!!$                   &set")')
!!$              write(*,'(" which is required for tetrahedron method")')
!!$              write(*,'(" commensurability tolerance: ",g18.10)') epscomm
!!$              write(*,'(" q-point (latt. coords.)   : ",3g18.10)') vql(:,iq)
!!$              write(*,'(" deviation                 : ",3g18.10)') vr/ngridk(:)
!!$              write(*,'(" minimum nonzero coords.   : ",3g18.10)')1.d0/ngridk(:)
!!$              write(*,*)
!!$              call terminate
!!$           end if
!!$           iqnr=1+iv(1)+ngridq(1)*iv(2)+ngridq(1)*ngridq(2)*iv(3)
           ! generate link array for tetrahedra
           call gentetlink(vql(1,iq))
        end if
        ! prefactor
        t1=-4.d0*pi/omega
        f(:,:)=t1*f(:,:)
        ! tetrahedron method
        forall (ik=1:nkpt,ist1=1:nstsv)
           e1(ist1,ik)=evalsv(ist1,ik)
        end forall
        ! scissors correction needed in input energies for tetrahedron method
        where (e1.gt.efermi) e1=e1+scissor
        do iw=1,nwdos
           if ((modulo(iw,nwdos/10).eq.0).or.(iw.eq.nwdos)) &
                write(*,'("Info(linopt): tetrahedron weights for ",I6," of ",&
                &I6," w-points")') iw,nwdos
           ! it seems that frequency should be non-zero for tetcw (?)
           ! see Ricardo's code
           if (abs(w(iw)).lt.epstetra) w(iw)=epstetra
           ! switch 2 below in tetcw defines bulk integration for real part
           call tetcw(nkpt,ntet,nstfv,wtet,e1,tnodes,link,tvol,efermi, &
                w(iw),2,cw)
           call tetcw(nkpt,ntet,nstfv,wtet,e1,tnodes,link,tvol,efermi, &
                -w(iw),2,cwa)
           ! switch 4 below in tetcw defines surface integration for imag. part
           call tetcw(nkpt,ntet,nstfv,wtet,e1,tnodes,link,tvol,efermi, &
                w(iw),4,cwsurf)
           ! summation using weights from tetrahedron method
           sum1=0.d0
           sum2=0.d0
           do ik=1,nkpt
              m=0
              do ist1=1,nstsv
                 do ist2=1,nstsv
                    m=m+1
                    t3=1.d0
                    if (.not.tqfmt) t3=1.d0/e(m,ik)**2
                    if ((tqfmt.and.intraband).or.(ist1.ne.ist2)) then
                       ! real part, resonant contribution
                       if (ist1.le.ist2) sum1=sum1+cw(ist1,ist2,ik)* &
                            f(m,ik)*t3
                       ! real part, anti-resonant contribution
                       if (ist1.gt.ist2) sum1=sum1-cwa(ist2,ist1,ik)* &
                            f(m,ik)*t3
                       ! imaginary part (only resonant contribution by theory
                       ! for positive frequencies)
                       if (ist1.le.ist2) sum2=sum2+cwsurf(ist1,ist2,ik)* &
                            f(m,ik)
                    end if
                 end do
              end do
           end do ! ik
           eps1r(iw)=sum1
           eps2(iw)=sum2
           ! divide by omega^2
           if (.not.tqfmt) then
              t1=w(iw)-scissor !SAG
              if (abs(t1).gt.eps) then
                 eps2(iw)=eps2(iw)/(t1**2)
              else
                 eps2(iw)=0.d0
              end if
           end if
        end do ! iw
#endif
     else if (optltz) then
        ! prefactor
        t1=-4.d0*pi/omega
        f(:,:)=t1*f(:,:)
        ! Lorentzian broadening
        do iw=1,nwdos
           sum1=0.d0
           sum2=0.d0
           do ik=1,nkpt
              m=0
              do ist1=1,nstsv
                 do ist2=1,nstsv
                    m=m+1
                    scis=0.d0
                    if ((evalsv(ist1,ik).le.efermi).and. &
                         (evalsv(ist2,ik).gt.efermi)) scis=-scissor
                    if ((evalsv(ist1,ik).gt.efermi).and. &
                         (evalsv(ist2,ik).le.efermi)) scis=scissor
                    t3=1.d0
                    if (.not.tqfmt) then
                       if (abs(e(m,ik)-scis).gt.eps) then
                          t3=1.d0/(e(m,ik)-scis)**2
                       else
                          t3=0.d0
                          if ((ist1.ne.ist2).and.(iw.eq.1)) then
                             write(*,'(a,3i8,g18.10)') 'divergent energy &
                                  &denominator for ik,ist1,ist2:',ik,ist1, &
                                  ist2,e(m,ik)
                          end if
                       end if
                    end if
                    if ((tqfmt.and.intraband).or.(ist1.ne.ist2)) then
                       sum1=sum1+wkpt(ik)* f(m,ik)* &
                            dble(1.d0/(e(m,ik)+w(iw)+zi*optswidth))*t3
                       sum2=sum2+wkpt(ik)* f(m,ik)* &
                            aimag(1.d0/(e(m,ik)+w(iw)+zi*optswidth))*t3
                    end if
                 end do
              end do
           end do ! ik
           eps1r(iw)=sum1
           eps2(iw)=sum2
        end do
     else ! if (tetra)
        ! prefactor
        t1=-4.d0*(pi**2)/omega
        f(:,:)=t1*f(:,:)
        ! calculate imaginary part of the interband dielectric function
        call brzint(nsmdos,ngridk,nsk,ikmap,nwdos,wdos,n,n,e,f,eps2)
        if (.not.tqfmt) then
           do iw=1,nwdos
              t1=w(iw)-scissor !SAG
              if (abs(t1).gt.eps) then
                 eps2(iw)=eps2(iw)/(t1**2)
              else
                 eps2(iw)=0.d0
              end if
           end do
        end if
        ! calculate the intraband Drude-like contribution and plasma frequency
        if (intraband.and.(.not.tqfmt)) then
           write(fname,'("PLASMA_",2I1,".OUT")') i1,i2
           open(62,file=trim(fname),action='WRITE',form='FORMATTED')
           wd(1)=efermi-0.001d0
           wd(2)=efermi+0.001d0
           call brzint(nsmdos,ngridk,nsk,ikmap,3,wd,nstsv,nstsv,evalsv,&
                pmatint,g)
           wplas=sqrt(g(2)*8.d0*pi/omega)
           do iw=1,nwdos
              if (abs(w(iw)).gt.eps) then
                 ! add the intraband contribution to the imaginary part of the
                 ! tensor
                 eps2(iw)=eps2(iw)+swidth*(wplas**2)/(w(iw)*(w(iw)**2+&
                      swidth**2))
              end if
           end do
           ! write plasma frequency to file
           write(62,'(G18.10," : plasma frequency")') wplas
           close(62)
        end if
#ifdef TETRA
     end if ! if (tetra)
#endif
     ! calculate real part of the dielectric function
     if (i1.eq.i2) then
        t1=1.d0
     else
        t1=0.d0
     end if
     ! Kramers-Kronig transform to find real part of dielectric tensor
     do iw=1,nwdos
        do jw=1,nwdos
           t2=w(jw)**2-w(iw)**2
           if (abs(t2).gt.eps) then
              fw(jw)=w(jw)*eps2(jw)/t2
           else
              fw(jw)=0.d0
           end if
        end do
        call fderiv(-1,nwdos,w,fw,g,cf)
        eps1(iw)=t1+(2.d0/pi)*g(nwdos)
        !<sag>
        if (bzsmpl.ge.1) eps1r(iw)=t1+eps1r(iw)
        !</sag>
     end do
     ! write dielectric function to a file
     !<sag action="modifications">
     escal=1.d0
     if (tev) escal=27.2114d0
     ! modified output variables and format
     do iw=1,nwdos
        select case(bzsmpl)
        case(0)
           write(60,'(3G18.10)') escal*w(iw),eps1(iw),eps2(iw)
        case(1,2)
           write(60,'(4G18.10)') escal*w(iw),eps1(iw),eps2(iw),eps1r(iw)
        end select
     end do
     !</sag>
     ! calculate optical conductivity
     sigma1(:)=(eps2(:))*w(:)/(4.d0*pi)
     sigma2(:)=-(eps1(:)-t1)*w(:)/(4.d0*pi)
     ! write the interband optical conductivity to a file
     do iw=1,nwdos
        write(61,'(2G18.10)') w(iw),sigma1(iw)
     end do
     write(61,'("     ")')
     do iw=1,nwdos
        write(61,'(2G18.10)') w(iw),sigma2(iw)
     end do
     close(60)
     close(61)
     ! end loop over number of components
  end do
  !<sag action="modifications">
  if (.not.tqfmt) close(50)
  !</sag>
  write(*,*)
  write(*,'("Info(linopt):")')
  write(*,'(" dielectric tensor written to EPSILON_ij.OUT")')
  write(*,*)
  write(*,'(" optical conductivity written to SIGMA_ij.OUT")')
  if (intraband) then
     write(*,*)
     write(*,'(" plasma frequency written to PLASMA_ij.OUT")')
  end if
  write(*,*)
  write(*,'(" for all requested components i,j")')
  deallocate(w,fw,g,cf,e,f,sc,d)
  deallocate(eps1,eps2)
  deallocate(sigma1,sigma2)
  deallocate(evecfv,evecsv,pmat,pmatint)
  if (usegdft) deallocate(delta,apwalm)
#ifdef TETRA
  if (tetra) then
     deallocate(e1,cw,cwa,cwsurf)
  end if
#endif
  if (bzsmpl.ge.1) then
     deallocate(eps1r)
  end if
  return
end subroutine linopt
