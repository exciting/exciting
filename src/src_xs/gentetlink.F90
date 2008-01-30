
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gentetlink
! !INTERFACE:
subroutine gentetlink(vpl)
! !USES:
  use modmain
  use modtetra
! !DESCRIPTION:
!   Generates an array connecting the tetrahedra of the $\mathbf{k}$-point with
!   the ones of the  $\mathbf{k}+\mathbf{q}$-point. Interface routine
!   referencing the {\tt libbzint} library of Ricardo Gomez-Abal.
!
! !REVISION HISTORY:
!   Created January 2008 (SAG)
!EOP
!BOC
  implicit none
  ! arguments
  real(8), intent(in) :: vpl(3)
  ! local variables
  real(8), parameter :: epscomm=1.d-5
  real(8) :: vr(3)
  integer :: j,iv(3),iqnr
  integer, allocatable :: ivkt(:,:),ivqt(:,:),tnodest(:,:),wtett(:)

!!$integer,allocatable::linkt(:,:)

  ! get index to reducible q-point which is commensurate to k-point set
  vr(:)=vpl(:)*ngridk(:)
  call r3frac(epslat,vr,iv)
  if (sum(abs(vr)).gt.epscomm) then
     write(*,*)
     write(*,'("Error(m_tetcalccwq): q-point not commensurate with k-point &
          &set")')
     write(*,'(" which is required for tetrahedron method")')
     write(*,'(" commensurability tolerance: ",g18.10)') epscomm
     write(*,'(" q-point (latt. coords.)   : ",3g18.10)') vpl
     write(*,'(" deviation                 : ",3g18.10)') vr/ngridk(:)
     write(*,'(" minimum nonzero coords.   : ",3g18.10)') 1.d0/ngridk(:)
     write(*,*)
     call terminate
  end if
  iqnr=1+iv(1)+ngridq(1)*iv(2)+ngridq(1)*ngridq(2)*iv(3)
  ! check if k-point set is not reduced for q-point different from Gamma point
  if ((nkpt.ne.nkptnr).and.(iqnr.ne.1)) then
     write(*,*)
     write(*,'("Error(gentetlink): k-point set is reduced by symmetries and &
          &q-point is not Gamma point")')
     write(*,*)
     call terminate
  end if
  ! allocate link array
  if (allocated(link)) deallocate(link)
  allocate(link(6*nkpt))
  ! quick return for Gamma q-point
  if (iqnr.eq.0) then
     forall (j=1:6*nkpt) link(j)=j
     return
  end if
  ! allocate local arrays
  allocate(ivkt(3,nkptnr),ivqt(3,nkptnr))
  allocate(wtett(6*nkptnr),tnodest(4,6*nkptnr))
  ! generate fraction for k-point offset
  call r3fraction(vkloff,ikloff,dkloff)



!SAG **************************************************************************
  ! call to libbzint-routine
  call kqgen_exciting(bvec,ngridk,ikloff,dkloff,nkpt,iqnr,ivkt,ivqt,dvk,dvq, &
       ntet,tnodes      ,wtet     ,link,tvol)


!!$allocate(linkt(6*nkptnr,nkptnr),kqid(nkpt,nkpt))
!!$  call kqgen(bvec,ngridk,ikloff,dkloff,nkpt, ivkt,ivqt,dvk,dvq,   &
!!$     &                 kqid,  ntet,tnodes     ,wtet       ,linkt,tvol)
!!$! assign link array
!!$link(:)=linkt(:,iqnr)
!!$deallocate(linkt,kqid)



  ! Note: we are using the weights and nodes from "kgen" called in "init1"
  ! deallocate local arrays
  deallocate(ivkt,ivqt)
  deallocate(wtett,tnodest)
end subroutine gentetlink
!EOC
