
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gentetlink
! !INTERFACE:
subroutine gentetlink(iqnr)
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
  integer, intent(in) :: iqnr
  ! local variables
  integer :: j
  integer, allocatable :: ivkt(:,:),ivqt(:,:),tnodest(:,:),wtett(:)

integer,allocatable::linkt(:,:)

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
!!$  ! call to libbzint-routine
!!$  call kqgen_exciting(bvec,ngridk,ikloff,dkloff,nkpt,iqnr,ivkt,ivqt,dvk,dvq, &
!!$       ntet,tnodest      ,wtett     ,link,tvol)


allocate(linkt(6*nkptnr,nkptnr),kqid(nkpt,nkpt))
  call kqgen(bvec,ngridk,ikloff,dkloff,nkpt, ivkt,ivqt,dvk,dvq,   &
     &                 kqid,  ntet,tnodes     ,wtet       ,linkt,tvol)

! assign link array
link(:)=linkt(:,iqnr)
deallocate(linkt,kqid)



  ! Note: we are using the weights and nodes from "kgen" called in "init1"
  ! deallocate local arrays
  deallocate(ivkt,ivqt)
  deallocate(wtett,tnodest)
end subroutine gentetlink
!EOC
