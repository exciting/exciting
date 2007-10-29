
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine init2xs
  use modmain
  use modtddft
  use modtetra
  use modmpi
  use m_genqvkloff
  use m_findkmapkq
  use m_getunit
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'init2xs'
  character(256) :: bname
  integer :: nsym,iq,l,m,lm
  logical :: existent

  !--------------------!
  !     file names     !
  !--------------------!
  call genfilname(basename='PMAT_TD',filnam=fnpmat)
  call genfilname(basename='PMAT_TD',procs=procs,rank=rank,filnam=fnpmat_t)

  !------------------------------------!
  !     angular momentum variables     !
  !------------------------------------!
  
  ! if not specified in input file set lmaxapwtd to lmaxmat
  if (lmaxapwtd == -1) lmaxapwtd=lmaxmat
  ! check lmaxapwtd
  if (lmaxapwtd > lmaxapw) then
     write(*,*) 'Error('//thisnam//'): lmaxapwtd > lmaxapw:', lmaxapwtd
     call terminate
  end if
  lmmaxapwtd = (lmaxapwtd+1)**2
  lmmaxemat = (lmaxemat+1)**2
  lmaxmax = maxval((/lmaxemat,lolmax, lmaxvr,lmaxapw,lmaxmat,lmaxinr/))
  lmmaxmax = maxval((/lmmaxemat,lolmmax, lmmaxvr,lmmaxapw,lmmaxmat,lmmaxinr/))

  ! index to (l,m) pairs (overall)
  if (allocated(idxxlm)) deallocate(idxxlm)
  allocate(idxxlm(0:lmaxmax,-lmaxmax:lmaxmax))
  idxxlm(:,:) = 0
  lm=0
  do l=0,lmaxmax
     do m=-l,l
        lm=lm+1
        idxxlm(l,m)=lm
     end do
  end do

  ! array of i**l values
  if (allocated(zmil)) deallocate(zmil)
  allocate(zmil(0:lmaxmax))
  do l=0,lmaxmax
     zmil(l)=(-zi)**l
  end do

  !---------------------!
  !     q-point set     !
  !---------------------!
  ! check input type of q-point
  if ((trim(qtype).eq.'grid').or.(trim(qtype).eq.'zero')) then
     if (trim(qtype).eq.'zero') then
        ngridq(:) = 1
        vqloff(:) = 0.d0
     end if
     if (task==400) then
        ngridq(:)=ngridk(:)
     end if
     if (allocated(ivq)) deallocate(ivq)
     allocate(ivq(3,ngridq(1)*ngridq(2)*ngridq(3)))
     if (allocated(vql)) deallocate(vql)
     allocate(vql(3,ngridq(1)*ngridq(2)*ngridq(3)))
     if (allocated(vqc)) deallocate(vqc)
     allocate(vqc(3,ngridq(1)*ngridq(2)*ngridq(3)))
     if (allocated(wqpt)) deallocate(wqpt)
     allocate(wqpt(ngridq(1)*ngridq(2)*ngridq(3)))
     if (allocated(iqmap)) deallocate(iqmap)
     allocate(iqmap(0:ngridq(1)-1,0:ngridq(2)-1,0:ngridq(3)-1))
     ! generate reduced q-point set
     call genppts(reduceq,ngridq,vqloff,nqpt,iqmap,ivq,vql,vqc,wqpt)
  else if (trim(qtype).eq.'list') then
     ! read q-points from list (file)
     inquire(file=trim(qlist),exist=existent)
     if (.not.existent) then
        write(*,*) 'Error('//thisnam//'): no qlist file: '//trim(qlist)
        call terminate
     end if
     call getunit(unit1)
     open(unit=unit1,file=trim(qlist),form='formatted',action='read', &
          status='old')
     ! read in the list with blockname "qlist"
10   continue
     read(unit1,*,end=20) bname
     ! check for a comment
     if ((scan(trim(bname),'!').eq.1).or.(scan(trim(bname),'#').eq.1)) goto 10
     select case(trim(bname))
     case('qlist')
        read(unit1,*) nqpt
        if (allocated(vql)) deallocate(vql)
        allocate(vql(3,nqpt))
        if (allocated(vqc)) deallocate(vqc)
        allocate(vqc(3,nqpt))
        if (allocated(wqpt)) deallocate(wqpt)
        allocate(wqpt(nqpt))
        wqpt(:) = 1.d0
        do iq=1,nqpt
           read(unit1,*) vql(:,iq)
           ! generate q-points in Cartesian coordinates
           vqc(:,iq)=vql(1,iq)*bvec(:,1)+vql(2,iq)*bvec(:,2)+ &
                vql(3,iq)*bvec(:,3)
        end do
     case('')
        goto 10
     case default
        write(*,*) 'Error('//thisnam//'): qlist-block not found: '//trim(bname)
        call terminate
     end select
     goto 10
20   continue
     close(unit1)
  else
     write(*,*) 'Error('//thisnam//'): unknown qtype: '//trim(qtype)
     call terminate
  end if

  ! check for Gamma point
  tq1gamma = .false.
  if (all(vql(:,1).eq.0)) tq1gamma = .true.

  ! find (little) group of q
  do iq=1,nqpt
!*******     call findlitgq(vql(:,iq),nsymcrysq,scqmap)
  end do

  !-----------------------!
  !     k+q-point set     !
  !-----------------------!
  if (allocated(qvkloff)) deallocate(qvkloff)
  if (allocated(ikmapikq)) deallocate(ikmapikq)
  allocate(qvkloff(3,nqpt))
  allocate(ikmapikq(nqpt,nkpt))
  do iq=1,nqpt
     ! offset for k+q-point set derived from q-point
     call genqvkloff(vql(:,iq),qvkloff(:,iq))
     ! map from k-point index to k+q point index for same k
     call findkmapkq(iq,vql(:,iq),qvkloff(:,iq),ikmapikq(iq,:))
  end do

  !---------------------!
  !     G+q-point set   !
  !---------------------!
  ! checking
  if (gqmax.ge.gkmax) then
     write(*,*) 'Warning('//thisnam//'): gqmax >= gkmax: '
  end if
  ! maximum number of G+q vectors for all q
  call getngqmax
  ! allocate the G+q-vector arrays
  if (allocated(ngq)) deallocate(ngq)
  allocate(ngq(nqpt))
  if (allocated(igqig)) deallocate(igqig)
  allocate(igqig(ngqmax,nqpt))
  if (allocated(vgql)) deallocate(vgql)
  allocate(vgql(3,ngqmax,nqpt))
  if (allocated(vgqc)) deallocate(vgqc)
  allocate(vgqc(3,ngqmax,nqpt))
  if (allocated(gqc)) deallocate(gqc)
  allocate(gqc(ngqmax,nqpt))
  if (allocated(tpgqc)) deallocate(tpgqc)
  allocate(tpgqc(2,ngqmax,nqpt))
  if (allocated(sfacgq)) deallocate(sfacgq)
  allocate(sfacgq(ngqmax,natmtot,nqpt))
  if (allocated(ylmgq)) deallocate(ylmgq)
  allocate(ylmgq(lmmaxvr,ngqmax,nqpt))
  do iq=1,nqpt
     ! generate G+q vectors
     call gengqvec(vql(1,iq),vqc(1,iq),ngq(iq),igqig(1,iq), &
          vgql(1,1,iq),vgqc(1,1,iq),gqc(1,iq),tpgqc(1,1,iq))
     ! generate structure factors for G-vectors
     call gensfacgp(ngq(iq),vgqc(1,1,iq),ngq(iq),sfacgq(1,1,iq))
     ! spherical harmonics for G+q-vectors
     call genylmgq(iq)
  end do

  !------------------------!
  !     radial functions   !
  !------------------------!
  ! read density and potentials from file
  call readstate
  ! find the new linearisation energies
  call linengy
  ! generate the APW radial functions
  call genapwfr
  ! generate the local-orbital radial functions
  call genlofr

  !--------------------------!
  !     occupation numbers   !
  !--------------------------!
  ! number of occupied valence states (valence band states)
  nstval = nstsv - nempty - 1
  ! number of unoccupied valence states (conduction band states)
  nstcon = nempty + 1

  !-------------------------------!
  !     analytic continuation     !
  !-------------------------------!
  ! if imaginary frequencies intervals are not specified
  if (nwacont.eq.0) nwacont = nwdos
  nwdf=nwdos
  if (acont) nwdf=nwacont

  !------------------------------------!
  !     response function variables    !
  !------------------------------------!
  if (allocated(mdfrpa)) deallocate(mdfrpa)
  allocate(mdfrpa(nwdos,3,2))
  mdfrpa(:,:,:)=0.d0
  if (allocated(mdfrpad)) deallocate(mdfrpad)
  allocate(mdfrpad(nwdos,3))
  mdfrpad(:,:)=0.d0

  !----------------------------!
  !     xc-kernel variables    !
  !----------------------------!
  if (allocated(fxc0)) deallocate(fxc0)
  allocate(fxc0(nwdos,3))
  fxc0(:,:)=0.d0
  if (allocated(fxc0d)) deallocate(fxc0d)
  allocate(fxc0d(nwdos,3))
  fxc0d(:,:)=0.d0

  !---------------------------!
  !     exciton variables     !
  !---------------------------!
  if (allocated(excite)) deallocate(excite)
  allocate(excite(nexcitmax,3))
  if (allocated(excito)) deallocate(excito)
  allocate(excito(nexcitmax,3))
  excite(:,:)=0.d0
  excito(:,:)=0.d0

  !---------------------------------!
  !     k-point parallelization     !
  !---------------------------------!
  if (procs > nkpt) then
     procs=nkpt
     write(*,*) 'Warning('//thisnam//'): procs > nkpt: resetting to nkpt'
     if (rank >= nkpt) then
        write(*,*) 'Warning('//thisnam//'): rank > nkpt: skipping this &
             &process - this should not happen within an MPI run'
        call tdepilog
        ! just stop current process
        stop
     end if
  end if

end subroutine init2xs
