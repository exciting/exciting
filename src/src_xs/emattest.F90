
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine emattest
  use modmain
  use modxs
  use m_getunit
  use m_getpmat
  use m_getemat
  use m_genfilname
  implicit none
  complex(8), allocatable :: pmat(:,:,:,:), x(:,:,:,:)
  real(8), allocatable :: d(:,:,:),scis12(:,:),scis21(:,:)
  ! compare first Q-point
  integer, parameter :: iq=1
  complex(8) :: x_sc, p_sc
  real(8) :: d1,d2,d3,a,p
  integer :: n,ik,ikq,ist1,ist2
  character(256) :: filename
  logical, external :: tqgamma
  call init0
  call init1
  call init2xs
  call xssave0
  n = ngq(iq)
  if (tqgamma(iq)) then
     write(*,*)
     write(*,'("Info(emattest): Q-point is Gamma point - no comparison &
          &possible")')
     write(*,*)
     call terminate
  end if
  ! file extension for q-point
  call genfilname(iqmt=iq,setfilext=.true.)
  ! calculate k+q and G+k+q related variables
  call init1xs(vql(1,iq))
  call findocclims(iq,istocc0,istocc,istunocc0,istunocc,isto0,isto,istu0,istu)
  ! allocate arrays
  if (allocated(deou)) deallocate(deou)
  if (allocated(deuo)) deallocate(deuo)
  allocate(deou(nst1,nst2))
  allocate(deuo(nst2,nst1))
  if (allocated(docc12)) deallocate(docc12)
  if (allocated(docc21)) deallocate(docc21)
  allocate(docc12(nst1,nst2))
  allocate(docc21(nst2,nst1))
  if (allocated(xiou)) deallocate(xiou)
  if (allocated(xiuo)) deallocate(xiuo)
  allocate(xiou(nst1,nst2,n))
  allocate(xiuo(nst2,nst1,n))
  ! allocate local arrays
  allocate(x(nst1,nst2,n,nkpt))
  allocate(d(nst1,nst2,nkpt))
  allocate(pmat(3,nstsv,nstsv,nkpt))
  allocate(scis12(nst1,nst2))
  allocate(scis21(nst2,nst1))
  call getunit(unit1)
  call genfilname(basename='emat_pmat',iqmt=iq,filnam=filename)
  open(unit1,file=trim(filename),action='write',status='replace')
  ! annotate magnitude of q-vector
  write(*,*) 'Info(emattest): length of q-vector:',gqc(1,iq)
  ! test matrix elements
  do ik=1,nkpt
     ! read matrix elemets of exponential expression
     call getpmat(ik,vkl0,1,nstsv,.true.,trim(fnpmat),pmat(:,:,:,ik))
     ! read matrix elemets of exponential expression
     call getemat(iq,ik,.true.,trim(fnemat),xiou,xiuo)
     ikq=ikmapikq(ik,iq)
     call getdevaldoccsv(iq,ik,ikq,istlo1,isthi1,istlo2,isthi2,deou,docc12, &
          scis12)
     call getdevaldoccsv(iq,ik,ikq,istlo2,isthi2,istlo1,isthi1,deuo,docc21, &
          scis21)
     x(:,:,:,ik) = xiou(:,:,:)
     d(:,:,ik) = deou
     do ist1=1,nst1
        do ist2=1,nst2
           x_sc = x(ist1,ist2,1,ik)/gqc(1,iq)
           p_sc = dot_product(vgqc(:,1,iq)/gqc(1,iq), &
                pmat(:,istlo1-1+ist1,istlo2-1+ist2,ik))/(-d(ist1,ist2,ik))
           a = dble(x_sc)**2 + aimag(x_sc)**2
           p = dble(p_sc)**2 + aimag(p_sc)**2
           d1 = abs(x_sc-p_sc)
           d2 = min(d1,abs(x_sc+p_sc))
           d3 = abs(x_sc) - abs(p_sc)
           write(unit1,'(100g18.10)') ik,ist1,ist2, &
                x_sc, p_sc, d1, d2, d3, a, p
        end do
     end do
  end do
  deallocate(deou,deuo,docc12,docc21,scis12,scis21,xiou,xiuo,pmat,x,d)
  close(unit1)
  close(unit3)
  close(unit4)
end subroutine emattest
