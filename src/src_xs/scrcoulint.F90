
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine scrcoulint
  use modmain
  use modmpi
  use modxs
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam='scrcoulint'
  integer :: iknr,jknr,iq,j
  real(8) :: vklofft(3),vq(3),v1(3),v2(3),s(3,3),t1
  integer :: ngridkt(3),iv(3),lspl
  logical :: nosymt,reducekt
  real(8), external :: r3taxi
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
  call init0
  call init1
  call init2xs
  emattype=2
  ! check number of empty states
  if (nemptyscr.lt.nempty) then
     write(*,*)
     write(*,'("Error(",a,"): too few empty states in screening eigenvector &
          &file - the screening should include many empty states &
          &(BSE/screening)",2i8)') trim(thisnam),nempty,nemptyscr
     write(*,*)
     call terminate
  end if
  call genfilname(dotext='_SCI.OUT',setfilext=.true.)
  if (rank.eq.0) call writekpts


  ! loop over non-reduced number of k-points
  do iknr=1,nkptnr     
     do jknr=iknr,nkptnr
        iv(:)=ivknr(:,jknr)-ivknr(:,iknr)
        iv(:)=modulo(iv(:),ngridk(:))
        ! q-point (non-reduced)
        iq=ikmap(iv(1),iv(2),iv(3))
        vq(:)=vql(:,iq)
        ! q-point (reduced)
        iqnr=ikmapnr(iv(1),iv(2),iv(3))
        vqnr(:)=vklnr(:,iqnr)-vkloff(:)
        ! symmetry that transforms non-reduced q-point to reduced one
        do isym=1,nsymcrys
           lspl=lsplsymc(isym)
           s(:,:)=dble(symlat(:,:,lspl))
           call r3mtv(s,v1,v2)
           call r3frac(epslat,v2,iv)
           t1=r3taxi(vklnr(1,iknr),v2)
           if (t1.lt.epslat) then ??????????

        end do

     end do
  end do





  ! restore global variables
  nosym=nosymt
  reducek=reducekt
  ngridk(:)=ngridkt(:)
  vkloff(:)=vklofft(:)
  write(unitout,'(a)') "Info("//trim(thisnam)//"): Screening finished"
end subroutine scrcoulint
