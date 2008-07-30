
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modtetra
  ! Variable names taken from the GW implementation into the EXCITING code
  ! version 0.9.52 by R. Gomez-Abal.
  implicit none

  !----------------------------!
  !     ordering variables     !
  !----------------------------!
  ! map from library k-point index to default k-point index
  integer, allocatable :: iktet2ik(:)
  ! reverse map
  integer, allocatable :: ik2iktet(:)

  !--------------------------------------!
  !     tetrahedron method variables     !
  !--------------------------------------!
  ! true if tetrahedron method is used
  logical :: tetra
  ! integer k-point offset
  integer(4) :: ikloff(3)
  ! k-point offset divisor
  integer(4) :: dkloff
  ! k-points common divisor
  integer(4) :: dvk
  ! Number of tetrahedra
  integer(4) :: ntet
  ! index of the k-points corresponding to the nodes of each tetrahedron
  integer(4), allocatable :: tnodes(:,:)
  ! weight of each tetrahedron.
  integer(4), allocatable :: wtet(:)
  ! volume of the tetrahedra relative to the BZ volume
  real(8) :: tvol
  ! parameter specifying smalles diagonal in generic tetrahedron
  integer(4) :: mnd

  !---------------------------------!
  !     q-dependent convolution     !
  !---------------------------------!
  ! number of the tetrahedra linked to by the corresponding q vector
  integer(4), allocatable :: link(:), kqid(:,:)
  ! q-points common divisor
  integer(4) dvq

contains

!BOP
! !ROUTINE: rtorat
! !INTERFACE:
  subroutine rtorat(eps,n,x,k,div)
! !DESCRIPTION: 
!   This subroutine factorizes the real coordinates of a vector {\bf x}.
!   The output is an integer vector {\bf k}, such that $k(i)/{\rm div}=x(i)$
!   and
!   $$ |x(i)-k(i)/{\rm div}| < {\rm eps} $$.
!
! !REVISION HISTORY:
!   Created July 2008 by Sagmeister
!EOP
!BOC
    implicit none
    ! arguments
    real(8), intent(in) :: eps
    integer(4), intent(in) :: n
    real(8), intent(in) :: x(n)
    integer(4), intent(out) :: div
    integer(4), intent(out) :: k(n)
    ! local variables
    integer :: maxint
    real(8) :: dx
    maxint=nint(1.d0/eps)
    do div=1,maxint
       k(:)=nint(dble(div)*x(:))
       dx=maxval(abs(dble(k)/dble(div)-x))
       if (dx.lt.eps) exit
    end do
    if (dx.ge.eps) then
       write(*,*)
       write(*,'("Error(modtetra:rtorat): factorization failed")')
       write(*,'(" maximum integer :",i12)') maxint
       write(*,'(" tolerance       :",g18.10)') eps
       write(*,'(" deviation       :",g18.10)') dx
       write(*,*)
       stop
    end if
    if (dx.gt.1.d-12) then
       write(*,*)
       write(*,'("Warning(modtetra:rtorat): small deviation in factorization")')
       write(*,'(" maximum deviation :",g18.10)') dx
       write(*,*)
    end if
  end subroutine rtorat
!EOC


  subroutine r3fraction(r,n,d)
    implicit none
    ! arguments
    real(8), intent(in) :: r(3)
    integer, intent(out) :: n(3),d
    ! parameters
    real(8), parameter :: eps=1.d-5,eps2=1.d-3
    call rtorat(eps,3,r,n,d)
    ! check factorization
    if ((sum(abs(r)).lt.eps2).and.(sum(abs(r)).gt.0.d0)) then
       write(*,*)
       write(*,'("Warning(modtetra:r3fraction): very small offset:")')
       write(*,'(" kgen and related routines might fail")')
       write(*,*)
    end if
  end subroutine r3fraction


  subroutine geniktetmap(eps,nppt,ngridp,vploff,vpllib,vpl,ipmap)
    implicit none
    ! arguments
    real(8), intent(in) :: eps
    integer, intent(in) :: nppt
    integer, intent(in) :: ngridp(3)
    real(8), intent(in) :: vploff(3)
    real(8), intent(in) :: vpl(3,nppt)
    real(8), intent(in) :: vpllib(3,nppt)
    integer, intent(in) :: ipmap(0:ngridp(1)-1,0:ngridp(2)-1,0:ngridp(3)-1)
    ! local variables
    integer :: ip,ipd,iv(3)
    if (allocated(iktet2ik)) deallocate(iktet2ik)
    allocate(iktet2ik(nppt))
    if (allocated(ik2iktet)) deallocate(ik2iktet)
    allocate(ik2iktet(nppt))
    do ip=1,nppt
       ! grid coordinates of library k-point
       iv(:)=nint(vpllib(:,ip)*ngridp-vploff(:))
       ! index in default p-point set
       ipd=ipmap(iv(1),iv(2),iv(3))
       ! map from library to default
       iktet2ik(ip)=ipd
       ! reverse map
       ik2iktet(ipd)=ip
       ! check maps
       if (sum(abs(vpl(:,ipd)-vpllib(:,ip))).ge.eps) then
          write(*,*)
          write(*,'("Error(modtetra:geniktetmap): k-point mapping between")')
          write(*,'(" set of library and default set failed:")')
          write(*,'(" library k-point       :",3g18.10)') vpllib(:,ip)
          write(*,'(" mapped default k-point:",3g18.10)') vpl(:,ipd)
          write(*,'(" map from library to default:",2i8)') ip,ipd
          write(*,*)
          stop
       end if
    end do
  end subroutine geniktetmap


  subroutine writeiktetmap(filext,nkpt)
    implicit none
    ! arguments
    character(*), intent(in) :: filext
    integer, intent(in) :: nkpt
    ! local variables
    integer :: ik,un
    un=771
    open(un,file='KTETMAP'//trim(filext),form='formatted',action='write', &
         status='replace')
    write(un,'(i9,a)') nkpt, ' : nkpt; k-point, iktet2ik below'
    do ik=1,nkpt
       write(un,'(2i9)') ik,iktet2ik(ik)
    end do
    close(un)
    open(un,file='KTETIMAP'//trim(filext),form='formatted',action='write', &
         status='replace')
    write(un,'(i9,a)') nkpt, ' : nkpt; k-point, ik2iktet below'
    do ik=1,nkpt
       write(un,'(2i9)') ik,ik2iktet(ik)
    end do
    close(un)
  end subroutine writeiktetmap


  subroutine genkpts_tet(filext,eps,bvec,maxsymcrys,nsymcrys,lsplsymc,symlat, &
       reducek,ngridk,vkloff,nkpt,ikmap,vkl,wkpt)
    implicit none
    ! arguments
    character(*), intent(in) :: filext
    real(8), intent(in) :: eps
    real(8), intent(in) :: bvec(3,3)
    integer, intent(in) :: maxsymcrys
    integer, intent(in) :: nsymcrys
    integer, intent(in) :: lsplsymc(maxsymcrys)
    integer, intent(in) :: symlat(3,3,48)
    logical, intent(in) :: reducek
    integer, intent(in) :: ngridk(3)
    real(8), intent(in) :: vkloff(3)
    integer, intent(in) :: nkpt
    integer, intent(in) :: ikmap(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1)
    real(8), intent(in) :: vkl(3,ngridk(1)*ngridk(2)*ngridk(3))
    real(8), intent(in) :: wkpt(ngridk(1)*ngridk(2)*ngridk(3))
    ! local variables
    integer :: isym,lspl,i1,i2,nsymcryst
    integer :: ik,ikd,nkptlib
    real(8) :: wkptlib
    integer, allocatable :: symc(:,:,:)
    integer, allocatable :: ivk(:,:)
    integer, allocatable :: indirkp(:)
    integer, allocatable :: iwkp(:)
    real(8), allocatable :: vkllib(:,:)
    if (.not.tetra) return
    if (nsymcrys.gt.48) then
       write(*,*)
       write(*,'("Error(modtetra:genkpts_tet): number of crystal symmetries > &
            &48:")')
       write(*,'(" This does not work with the k-point generation of")')
       write(*,'(" the linear tetrahedron method.")')
       write(*,*)
       stop
    end if
    ! switch to exciting interface
    !call tetrasetifc('exciting')
    ! suppress debug output in tetrahedron integration library (0)
    call tetrasetdbglv(0)
    ! safer pointer handling in tetrahedron integration library (1)
    call tetrasetpointerhandling(1)
    ! set resonance type (1...resonant weights)
    call tetrasetresptype(1)
    ! set treatment of q-shifted k-mesh
    call tetrasetkplusq(.true.)
    ! report interface parameters
    call tetrareportsettings
    ! generate fraction for k-point offset
    call rtorat(1.d-4,3,vkloff,ikloff,dkloff)
    ! get rotational part of crystal symmetries 
    allocate(symc(3,3,nsymcrys))
    do isym=1,nsymcrys
       lspl=lsplsymc(isym)
       ! transpose of rotation for use with the library
       do i1=1,3
          do i2=1,3
             symc(i1,i2,isym)=symlat(i2,i1,lspl)
          end do
       end do
    end do    
    nsymcryst=1
    if (reducek) nsymcryst=nsymcrys
    ! allocate local variables
    allocate(indirkp(ngridk(1)*ngridk(2)*ngridk(3)))
    allocate(iwkp(ngridk(1)*ngridk(2)*ngridk(3)))
    allocate(ivk(3,ngridk(1)*ngridk(2)*ngridk(3)))
    allocate(vkllib(3,ngridk(1)*ngridk(2)*ngridk(3)))
    ! allocate weights
    if (allocated(wtet)) deallocate(wtet)
    allocate(wtet(1:ngridk(1)*ngridk(2)*ngridk(3)*6))
    wtet(:)=0
    ! allocate nodes
    if (allocated(tnodes)) deallocate(tnodes)
    allocate(tnodes(1:4,1:ngridk(1)*ngridk(2)*ngridk(3)*6))
    tnodes(:,:)=0
    ! generate k-point set by using library-routine
    call kgen(bvec,nsymcryst,symc,ngridk,ikloff,dkloff,nkptlib,ivk,dvk,indirkp,&
         iwkp,ntet,tnodes,wtet,tvol,mnd)
    if (nkptlib.ne.nkpt) then
       write(*,*)
       write(*,'("Error(modtetra:genkpts_tet): k-point set inconsistency for &
            &tetrahedron method")')
       write(*,'(" differring number of k-points (library/default)",2i8)') &
            nkptlib,nkpt
       write(*,*)
       stop
    end if
    ! k-point in lattice coordinates
    do ik=1,nkpt
       vkllib(:,ik)=dble(ivk(:,ik))/dble(dvk)
    end do
    ! generate map between the k-point set of the library and the default one
    call geniktetmap(eps,nkpt,ngridk,vkloff,vkllib,vkl,ikmap)
    ! check weights of k-points
    do ik=1,nkpt
       ikd=iktet2ik(ik)
       wkptlib=iwkp(ik)/dble(ngridk(1)*ngridk(2)*ngridk(3))
       if (abs(wkpt(ikd)-wkptlib).gt.eps) then
          write(*,*)
          write(*,'("Error(modtetra:genkpts_tet): differring weights:")')
          write(*,'(" k-point (default)",i8,3g18.10)') ikd,vkl(:,ikd)
          write(*,'(" weight (default)",g18.10)') wkpt(ikd)
          write(*,'(" weight (library)",g18.10)') wkptlib
          write(*,*)
          stop 
       end if
    end do
    deallocate(indirkp,iwkp,ivk,vkllib)
    ! write maps
    call writeiktetmap(filext,nkpt)
  end subroutine genkpts_tet


!BOP
! !ROUTINE: gentetlink
! !INTERFACE:
  subroutine gentetlink(vpl,tqw,eps,bvec,ngridk,vkloff,nkpt,nkptnr,vklnr, &
       ikmapnr)
! !DESCRIPTION:
!   Generates an array connecting the tetrahedra of the $\mathbf{k}$-point with
!   the ones of the  $\mathbf{k}+\mathbf{q}$-point. Interface routine
!   referencing the {\tt libbzint} library of Ricardo Gomez-Abal.
!
! !REVISION HISTORY:
!   Created January 2008 (Sagmeister)
!EOP
!BOC
    implicit none
    ! arguments
    real(8), intent(in) :: vpl(3)
    integer, intent(in) :: tqw
    real(8), intent(in) :: eps
    real(8), intent(in) :: bvec(3,3)
    integer, intent(in) :: ngridk(3)
    real(8), intent(in) :: vkloff(3)
    integer, intent(in) :: nkpt
    integer, intent(in) :: nkptnr
    real(8), intent(in) :: vklnr(3,ngridk(1)*ngridk(2)*ngridk(3))
    integer, intent(in) :: ikmapnr(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1)
    ! local variables
    real(8), parameter :: epscomm=1.d-5
    real(8) :: vr(3)
    integer :: j,iv(3),iqnr
    logical :: tqg
    integer, allocatable :: ivkt(:,:),ivqt(:,:),tnodest(:,:),wtett(:)
    real(8), external :: r3taxi
    tqg=sum(abs(vpl)).lt.eps
    ! get index to reducible q-point which is commensurate to k-point set
    vr(:)=vpl(:)*ngridk(:)
    call r3frac(eps,vr,iv)
    if (sum(abs(vr/ngridk)).gt.epscomm) then
       write(*,*)
       write(*,'("Error(gentetlink): q-point not commensurate with k-point &
            &set")')
       write(*,'(" which is required for tetrahedron method")')
       write(*,'(" commensurability tolerance: ",g18.10)') epscomm
       write(*,'(" q-point (latt. coords.)   : ",3g18.10)') vpl
       write(*,'(" deviation                 : ",3g18.10)') vr/ngridk(:)
       write(*,'(" minimum nonzero coords.   : ",3g18.10)') 1.d0/ngridk(:)
       write(*,*)
       call terminate
    end if
    iqnr=ikmapnr(iv(1),iv(2),iv(3))
    ! cross check q-point again
    vr(:)=vklnr(:,iqnr)-vkloff(:)/ngridk(:)
    if (abs(r3taxi(vpl,vr)).gt.epscomm) then
       write(*,*)
       write(*,'("Error(gentetlink): specified q-point does not match derived &
            &q-point on grid")')
       write(*,'(" specified q-point :",3g18.10)') vpl
       write(*,'(" derived q-point   :",3g18.10)') vr
       write(*,'(" non-reduced index :",i6)') iqnr
       write(*,*)
       call terminate
    end if
    write(*,*)
    write(*,'("Info(gentetlink): q-point on grid")')
    write(*,'(" q-point           :",3g18.10)') vr
    write(*,'(" non-reduced index :",i6)') iqnr
    write(*,*)
    ! check if k-point set is not reduced for q-point different from Gamma point
    if ((nkpt.ne.nkptnr).and.(.not.tqg)) then
       write(*,*)
       write(*,'("Error(gentetlink): k-point set is reduced by symmetries and &
            &q-point is not Gamma point")')
       write(*,*)
       call terminate
    end if
    ! allocate link array
    if (allocated(link)) deallocate(link)
    allocate(link(6*nkptnr))
    ! quick return for Gamma q-point
    if (tqg) then
       forall (j=1:6*nkptnr) link(j)=j
       return
    end if
    ! allocate local arrays
    allocate(ivkt(3,nkptnr),ivqt(3,nkptnr))
    allocate(wtett(6*nkptnr),tnodest(4,6*nkptnr))
    ! generate fraction for k-point offset
    call r3fraction(vkloff,ikloff,dkloff)
    ! call to library routine (generate link array, nodes and weights)
    call kqgen_exciting(bvec,ngridk,ikloff,dkloff,nkpt,iqnr,ivkt,ivqt,dvk,dvq, &
         ntet,tnodest,wtett,link,tvol)
    if (tqw.ne.0) then
       ! use weights and nodes from kqgen-routine
       tnodes(:,:)=tnodest(:,:)
       wtet(:)=wtett(:)
    else if (sum(abs(vpl)).gt.eps) then
       write(*,*)
       write(*,'("Warning(gentetlink): using WTET and TNODES arrays from KGEN &
            &routine")')
       write(*,'(" but arrays from KQGEN_EXCITING are differring for")')
       write(*,'(" non-Gamma q-point: ",3g18.10)') vpl
       write(*,*)
    end if
    ! deallocate local arrays
    deallocate(ivkt,ivqt)
    deallocate(wtett,tnodest)
  end subroutine gentetlink
!EOC


  subroutine fermitetifc(nkpt,nst,eval,chgval,spinpol,efermi,fermidos)
    implicit none
    ! arguments
    integer, intent(in) :: nkpt
    integer, intent(in) :: nst
    real(8), intent(in) :: eval(nst,nkpt)
    real(8), intent(in) :: chgval
    logical, intent(in) :: spinpol
    real(8), intent(out) :: efermi
    real(8), intent(out) :: fermidos
    ! local variables
    integer :: ik,ikd
    real(8), allocatable :: evallib(:,:)
    allocate(evallib(nst,nkpt))
    ! reorder energies to library order
    do ik=1,nkpt
       evallib(:,ik)=eval(:,iktet2ik(ik))
    end do
    ! call to library routine
    call fermitet(nkpt,nst,evallib,ntet,tnodes,wtet,tvol,chgval,spinpol, &
         efermi,fermidos,.false.)
    deallocate(evallib)
  end subroutine fermitetifc


  subroutine tetiwifc(nkpt,nst,eval,efermi,occ)
    implicit none
    ! arguments
    integer, intent(in) :: nkpt
    integer, intent(in) :: nst
    real(8), intent(in) :: eval(nst,nkpt)
    real(8), intent(in) :: efermi
    real(8), intent(out) :: occ(nst,nkpt)
    ! local variables
    integer :: ik,ikd
    real(8), allocatable :: evallib(:,:),occt(:)
    allocate(evallib(nst,nkpt),occt(nst))
    ! reorder energies to library order
    do ik=1,nkpt
       evallib(:,ik)=eval(:,iktet2ik(ik))
    end do
    ! call to library routine
    call tetiw(nkpt,ntet,nst,evallib,tnodes,wtet,tvol,efermi,occ)
    ! reorder occupation numbers to default order
    do ik=1,nkpt
       ikd=iktet2ik(ik)
       occt(:)=occ(:,ikd)
       occ(:,ikd)=occ(:,ik)
       occ(:,ik)=occt(:)
    end do
    deallocate(evallib,occt)
  end subroutine tetiwifc


  subroutine tetcwifc(nkpt,nst,eval,efermi,w,ifreq,cw)
    implicit none
    ! arguments
    integer, intent(in) :: nkpt
    integer, intent(in) :: nst
    real(8), intent(in) :: eval(nst,nkpt)
    real(8), intent(in) :: efermi
    real(8), intent(in) :: w
    integer, intent(in) :: ifreq       
    real(8), intent(out) :: cw(nst,nst,nkpt)
    ! local variables
    integer :: ik,ikd
    real(8), allocatable :: evallib(:,:),cwt(:,:)
    allocate(evallib(nst,nkpt),cwt(nst,nst))
    ! reorder energies to library order
    do ik=1,nkpt
       evallib(:,ik)=eval(:,iktet2ik(ik))
    end do
    ! call to library routine
    call tetcw(nkpt,ntet,nst,wtet,evallib,tnodes,link,tvol,efermi,w,ifreq,cw)
    ! reorder convolution weights to default order
    do ik=1,nkpt
       ikd=iktet2ik(ik)
       cwt(:,:)=cw(:,:,ikd)
       cw(:,:,ikd)=cw(:,:,ik)
       cw(:,:,ik)=cwt(:,:)
    end do
    deallocate(evallib,cwt)
  end subroutine tetcwifc

end module modtetra
