
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine geniscreen(iqr,nmax,n,scri)
  use modmain
  use modxs
  use invert
  use m_getunit
  implicit none
  ! arguments
  integer, intent(in) :: iqr,n,nmax
  complex(8), intent(out) :: scri(nmax,nmax)
  ! local variables
  complex(8), allocatable :: tm(:,:),tmi(:,:),b(:,:),bi(:,:),u(:,:),s(:,:)
  complex(8), allocatable :: scrn(:,:),scrnw(:,:,:),scrnh(:), s3(:,:),s3i(:,:)
  complex(8) :: scrnh0(3),scrnih0(3),f(3,3),l(3,3)
  integer :: oct,oct2,j,ig1,ig2,un
  integer, external :: octmap

  allocate(scrn(n,n),scrnw(n,2,3),scrnh(9),tm(n,n),tmi(n,n),s3(n+2,n+2), &
       s3i(n+2,n+2))
  ! read screening from file
  call getscreen(iqr,n,scrnh,scrnw,scrn)

  ! invert dielectric matrix for q-point going to zero in x-, y-, and
  ! z-direction for q=0
  if (iqr.eq.1) then
     ! head
     do oct=1,3
        do oct2=1,3
           j=octmap(oct,oct2)
           f(oct,oct2)=scrnh(j)
        end do
     end do
     if (n.gt.1) then
        allocate(b(n-1,n-1),bi(n-1,n-1),u(n-1,3),s(n-1,3))
        ! body
        b(:,:)=scrn(2:,2:)
        ! wings
        u(:,:)=scrnw(2:,1,:)
        ! invert body
        call zinvert_hermitian(scrherm,b,bi)
        s=matmul(bi,u)
        l=f-matmul(conjg(transpose(u)),s)
        deallocate(b,bi,u,s)
     else
        l=f
     end if

     ! write dielectric tensor to file
     call getunit(un)
     open(un,file='DIELTENSOR'//trim(filext),form='formatted',action='write', &
          status='replace')
     write(un,*)
     write(un,'(" dielectric tensor, real part:")')
     write(un,'(3f12.8)') (dble(l(oct,:)),oct=1,3)
     write(un,*)
     write(un,'(" dielectric tensor, imaginary part:")')
     write(un,'(3f12.8)') (aimag(l(oct,:)),oct=1,3)
     close(un)
     
     ! average the angular dependence of the inverse dielectric matrix at the
     ! Gamma q-point
     call angavdm(l,s,bi,scri)


  else
     ! case q/=0
     tm(:,:)=scrn(:,:)

     write(*,'(a,i5,2g18.10)') 'iq,    1/eps_00(q)             :', &
          iqr,1.d0/scrn(1,1)
     write(*,'(a,i5,2g18.10)') 'iq,      eps_00(q)             :', &
          
          iqr,scrn(1,1)
     call zinvert_hermitian(scrherm,tm,tmi)
     scri(1:n,1:n)=tmi(:,:)
  end if

  write(*,'(a,i5,2g18.10)') 'diel. matr.: iq,   eps^-1_00(q):', &
       iqr,scri(1,1)
  write(*,'(a,i5,2g18.10)') 'diel. matr.: iq, 1/eps^-1_00(q):', &
       iqr,1.d0/scri(1,1)
  write(*,*)

  deallocate(scrn,scrnw,scrnh,tm,tmi,   s3,s3i)
end subroutine geniscreen



subroutine angavdm(n,l,s,bi,av)
  use modmain
  use modxs
  implicit none
  ! arguments
  integer, intent(in) :: avtype
  integer, intent(in) :: n
  complex(8), intent(in) :: bi(n-1,n-1),s(n-1,3),l(3,3)
  complex(8), intent(out) :: av(n,n)
  ! local variables
  real(8), allocatable :: a(:,:),plat(:,:),p(:),m00(:),m0x(:),mxx(:)
  real(8), allocatable :: ei00(:),ei0x(:,:)
  allocate(a(3,lmmaxvr),plat(3,lmmaxvr),p(lmmaxvr))
  allocate(m00(lmmaxvr),m0x(lmmaxvr),mxx(lmmaxvr))
  allocate(ei00(lmmaxvr),ei0x(lmmaxvr,n-1),ei00lm(lmmaxvr),ei0xlm(lmmaxvr,n-1))
  a(:,:)=sphcov(:,:)
  ! unit vectors of spherical covering set in lattice coordinates
  plat=matmul(binv,a)
  ! distances to unit cell boundary in reciprocal space
  p(:)=1.d0/(2.d0*maxval(abs(plat),1))
  ! head which is the 1/(p*L*p) factor
  ei00(:)=1.d0/matmul(transpose(a),matmul(l,a))
  ! wings: -p*s * ei00
  do i=1,n-1
     ei0x(:,i)=-ei00(:)*matmul(s,a)
  end do

  ! forward transform shape dependent factor for head to spherical harmonics
  m00(:)=matmul(zfshtvr,p)
  ! forward transform shape dependent factor for wings to spherical harmonics
  m0x(:)=matmul(zfshtvr,(1.d0/2.d0)*p**2)
  ! forward transform shape dependent factor for body to spherical harmonics
  mxx(:)=matmul(zfshtvr,(1.d0/3.d0)*p**3)
  
  !forward transform to spherical harmonics
  ei00lm(:)=matmul(zfshtvr,ei00)
  ei0xlm(:)=matmul(zfshtvr,ei0x)

  ! angular average
  av(1,1)=sum(ei00lm*conjg(m00))
  

  deallocate(e,plat,p,m00,m0x,mxx,ei00,ei0x)
end subroutine angavdm





















