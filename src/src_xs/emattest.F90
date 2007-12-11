
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine emattest
  use modmain
  use modxs
  use m_getunit
  use m_getpmat
  use m_getemat
  use m_getdevalsv
  use m_genfilname
  implicit none
  complex(8), allocatable :: pmat(:,:,:,:), x(:,:,:,:)
  real(8), allocatable :: d(:,:,:)
  complex(8) :: x_sc, p_sc
  real(8) :: denom, d1,d2,d3,a,p
  integer :: iq,ik,istv,istc,n
  integer :: recl
  character(256) :: filename

  call init0
  call init1
  call init2xs
  call tdsave0

  ! SECOND Q-POINT not equal to zero
  iq=2
  ! file extension for q-point
  call genfilname(iq=iq,setfilext=.true.)
  ! shift k-mesh by q-point
  vkloff(:) = vkloff(:) + vql(:,iq)*ngridk(:)

  ! calculate k+q and G+k+q related variables
  call init1td

  n = ngq(iq)

  ! allocate arrays for eigenvalue differences
  if(allocated(deou)) deallocate(deou)
  if(allocated(deuo)) deallocate(deuo)
  allocate(deou(nstval,nstcon))
  allocate(deuo(nstcon,nstval))
  allocate(d(nstval,nstcon,nkpt))

  ! allocate matrix elements array
  if (allocated(xiou)) deallocate(xiou)
  if (allocated(xiuo)) deallocate(xiuo)
  allocate(xiou(nstval,nstcon,n))
  allocate(x(nstval,nstcon,n,nkpt))
  allocate(xiuo(nstcon,nstval,n))
  allocate(pmat(3,nstsv,nstsv,nkpt))


  call getunit(unit1)
  call genfilname(basename='emat_pmat',iq=iq,filnam=filename)
  open(unit1,file=trim(filename),action='write',status='replace')

  ! annotate magnitude of q-vector
  write(*,*) 'Info(emattest): length of q-vector:',gqc(1,iq)

  ! test matrix elements
  do ik=1,nkpt
     ! read matrix elemets of exponential expression
     call getpmat(ik,vkl0,.true.,trim(fnpmat),pmat(:,:,:,ik))
     ! read matrix elemets of exponential expression
     call getemat(iq,ik,.true.,trim(fnemat),xiou,xiuo)
     ! read Kohn-Sham eigenvalue differences
     call getdevalsv(iq,ik,.true.,trim(fndevalsv),deou,deuo)
     x(:,:,:,ik) = xiou(:,:,:)
     d(:,:,ik) = deou
     do istv=1,nstval
        do istc=1,nstcon
           x_sc = x(istv,istc,1,ik)/gqc(1,iq)
           p_sc = dot_product(vgqc(:,1,iq)/gqc(1,iq), &
                pmat(:,istv,nstval+istc,ik))/(-d(istv,istc,ik))
           a = dble(x_sc)**2 + aimag(x_sc)**2
           p = dble(p_sc)**2 + aimag(p_sc)**2
           d1 = abs(x_sc-p_sc)
           d2 = min(d1,abs(x_sc+p_sc))
           d3 = abs(x_sc) - abs(p_sc)
           write(unit1,'(100g18.10)') ik,istv,istc, &
                x_sc, p_sc, d1, d2, d3, a, p
        end do
     end do
  end do

  close(unit1)
  close(unit3)
  close(unit4)

end subroutine emattest
