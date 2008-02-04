
! Copyright (C) 2006-2007 S. Sagmeister, J. K. Dewhurst, S. Sharma and 
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine init2xs
  use modmain
  use modmpi
  use modxs
  use modtetra
  use m_getunit
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'init2xs'
  real(8) :: v(3)
  integer :: iq,l,m,lm,iv(3)

  !------------------------------------!
  !     angular momentum variables     !
  !------------------------------------!
  ! if not specified in input file set lmaxapwtd to lmaxmat
  if (lmaxapwtd.eq.-1) lmaxapwtd=lmaxmat
  ! check lmaxapwtd
  if (lmaxapwtd.gt.lmaxapw) then
     write(*,*) 'Error('//thisnam//'): lmaxapwtd > lmaxapw:', lmaxapwtd
     call terminate
  end if
  ! check lmaxemat
  if (lmaxemat.gt.lmaxapw) then
     write(*,*) 'Error('//thisnam//'): lmaxemat > lmaxapw:', lmaxemat
     call terminate
  end if
  if (lmaxemat.gt.lmaxapwtd) then
     write(*,*) 'Warning('//thisnam//'): lmaxemat > lmaxapwtd:', lmaxemat
     call terminate
  end if
  lmmaxapwtd=(lmaxapwtd+1)**2
  lmmaxemat=(lmaxemat+1)**2

  !---------------------!
  !     q-point set     !
  !---------------------!
  ! Q-/q-point set should have no offset
  vqloff(:)=0.d0
  ! assign momentum transfer Q-point set to q-point set
  if ((task.ge.300).and.(task.le.399)) then
     nqpt=nqptmt
     if (allocated(vqlmt)) deallocate(vqlmt)
     allocate(vqlmt(3,nqpt))
     if (allocated(ivgmt)) deallocate(ivgmt)
     allocate(ivgmt(3,nqpt))
     if (allocated(vql)) deallocate(vql)
     allocate(vql(3,nqpt))
     if (allocated(vqc)) deallocate(vqc)
     allocate(vqc(3,nqpt))
     do iq=1,nqpt
        v(:)=vgqlmt(:,iq)
        iv(:)=0
        ! map Q-point to Brillouin zone
        if (mdfqtype.eq.1) call r3frac(epslat,v,iv)
        vqlmt(:,iq)=v(:)
        ivgmt(:,iq)=iv(:)
        vql(:,iq)=vqlmt(:,iq)
        vqc(:,iq)=vql(1,iq)*bvec(:,1)+vql(2,iq)*bvec(:,2)+ &
             vql(3,iq)*bvec(:,3)
     end do
  end if
  ! generate q-point set from grid
  if ((task.ge.400).and.(task.le.499)) then
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
  end if

  ! find (little/small) group of q
  if (allocated(nsymcrysq)) deallocate(nsymcrysq)
  if (allocated(scqmap)) deallocate(scqmap)
  if (allocated(ivscwrapq)) deallocate(ivscwrapq)
  allocate(nsymcrysq(nqpt))
  allocate(scqmap(nsymcrys,nqpt))
  allocate(ivscwrapq(3,nsymcrys,nqpt))
  ! debug output
  if (dbglev.gt.1) then
     write(*,'(a)') 'Debug(init2xs):'
  end if
  do iq=1,nqpt
     call findgroupq(vql(1,iq),epslat,symlat,nsymcrys,lsplsymc,&
          nsymcrysq(iq),scqmap(1,iq),ivscwrapq(1,1,iq))
     ! debug output
     if (dbglev.gt.0) then
        write(*,'(a,i6,3g18.10,i6)') ' iq,vql,nsymcrysq:',iq,vql(:,iq),&
             nsymcrysq(iq)
        write(*,'(a)') ' scqmap,ivscwrapq below:'
        do l=1,nsymcrysq(iq)
           write(*,'(2i6,3x,3i6)') l,scqmap(l,iq),ivscwrapq(:,l,iq)
        end do
        write(*,*)
     end if
  end do

  !-----------------------!
  !     k+q-point set     !
  !-----------------------!
  if (allocated(qvkloff)) deallocate(qvkloff)
  allocate(qvkloff(3,0:nqpt))
  if (allocated(ikmapikq)) deallocate(ikmapikq)
  allocate(ikmapikq(nkpt,nqpt))
  qvkloff(:,0)=vkloff(:)
  do iq=1,nqpt
     ! offset for k+q-point set derived from q-point
     call genqvkloff(vql(1,iq),qvkloff(1,iq))
     ! map from k-point index to k+q point index for same k
     call findkmapkq(iq,vql(1,iq),qvkloff(1,iq),ikmapikq(1,iq))
  end do

  !---------------------!
  !     G+q-point set   !
  !---------------------!
  ! checking
  if (gqmax.ge.gkmax) then
     write(*,'(a,2g18.10)') 'Warning('//thisnam//'): gqmax >= gkmax: ',gqmax, &
          gkmax
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
  allocate(ylmgq(lmmaxapw,ngqmax,nqpt))
  if (allocated(ivgigq)) deallocate(ivgigq)
  allocate(ivgigq(intgqv(1,1):intgqv(1,2),intgqv(2,1):intgqv(2,2), &
       intgqv(3,1):intgqv(3,2),nqpt))
  do iq=1,nqpt
     ! generate G+q vectors
     call gengqvec(iq,vql(1,iq),vqc(1,iq),ngq(iq),igqig(1,iq), &
          vgql(1,1,iq),vgqc(1,1,iq),gqc(1,iq),tpgqc(1,1,iq))
     ! generate structure factors for G-vectors
     call gensfacgp(ngq(iq),vgqc(1,1,iq),ngq(iq),sfacgq(1,1,iq))
     ! spherical harmonics for G+q-vectors
     call genylmgq(iq,lmaxvr)
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
  if (allocated(occsv0)) deallocate(occsv0)
  allocate(occsv0(nstsv,nkpt))
  if (allocated(isto0)) deallocate(isto0)
  allocate(isto0(nkpt))
  if (allocated(isto)) deallocate(isto)
  allocate(isto(nkpt))
  if (allocated(istu0)) deallocate(istu0)
  allocate(istu0(nkpt))
  if (allocated(istu)) deallocate(istu)
  allocate(istu(nkpt))

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
  if (procs.gt.nkpt) then
     procs=nkpt
     write(*,*) 'Warning('//thisnam//'): procs > nkpt: resetting to nkpt'
     if (rank.ge.nkpt) then
        write(*,*) 'Warning('//thisnam//'): rank > nkpt: skipping this &
             &process - this should not happen within an MPI run'
        call xsfinit
        ! just stop current process
        stop
     end if
  end if

end subroutine init2xs
