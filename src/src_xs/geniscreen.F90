
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
  complex(8), allocatable :: scrn(:,:),scrnw(:,:,:),scrnh(:,:), s3(:,:),s3i(:,:)
  complex(8) :: scrnh0(3),scrnih0(3),f(3,3),l(3,3)
  integer :: un,oct,oct2,i,j,ig1,ig2

  allocate(scrn(n,n),scrnw(n,2,3),scrnh(3,3),tm(n,n),tmi(n,n),s3(n+2,n+2),s3i(n+2,n+2))
  ! read screening from file
  call getscreen(iqr,n,scrnh,scrnw,scrn)

  ! invert dielectric matrix for q-point going to zero in x-, y-, and
  ! z-direction for q=0
  if (iqr.eq.1) then
     ! head
     do oct=1,3
        do oct2=1,3
           f(oct,oct2)=scrnh(oct,oct2)
        end do
     end do
     allocate(b(n-1,n-1),bi(n-1,n-1),u(n-1,3),s(n-1,3))
     if (n.gt.1) then
        ! body
        b(:,:)=scrn(2:,2:)
        ! wings
        u(:,:)=conjg(scrnw(2:,1,:))
        ! invert body (optionally including Hermitian average)
        call zinvert_hermitian(scrherm,b,bi)
        s=matmul(bi,u)
        l=f-matmul(conjg(transpose(u)),s)

  write(*,*) 'f',f
  write(*,*)
  write(*,*) 'l',l
  write(*,*)

     else
        l=f
     end if

     ! write dielectric tensor to file
     call getunit(un)
     open(un,file='DIELTENSOR'//trim(filext),form='formatted',action='write', &
          status='replace')
     write(un,*)
     write(un,'(" dielectric tensor, real part:")')
     write(un,'(3f14.8)') (dble(l(oct,:)),oct=1,3)
     write(un,*)
     write(un,'(" dielectric tensor, imaginary part:")')
     write(un,'(3f14.8)') (aimag(l(oct,:)),oct=1,3)
     close(un)

     do oct=1,3
        tm(1,1)=f(oct,oct)
        tm(1,2:)=conjg(u(:,oct))
        tm(2:,2:)=b(:,:)
        ! set lower triangle
        do i=1,n
           do j=i+1,n
              tm(j,i)=conjg(tm(i,j))
           end do
        end do
        call zinvert_lapack(tm,tmi)
        write(7700+oct,'(2i8,3g18.10)') &
	 ((i,j,tmi(i,j),abs(tmi(i,j))**2,j=1,n),i=1,n)
     end do
!***************************************************************************


     do oct=1,3
        do oct2=1,3
           s3(oct,oct2)=scrnh(oct,oct2)
        end do
        if (n.gt.1) then
           s3(oct,4:)=scrnw(2:,1,oct)
        end if
   
        scrn(1,1)=scrnh(oct,oct)
        ! keep head of dielectric matrix for q=0
        scrnh0(oct)=scrn(1,1)
        if (n.gt.1) then
           scrn(1,2:n)=scrnw(2:n,1,oct)
           scrn(2:n,1)=scrnw(2:n,2,oct)
        end if
        ! store screening in temporary array
        tm=scrn

write(*,'(a,i5,2g18.10)') 'optcomp,    1/eps_00(q=0)      :', &
             oct,1.d0/scrnh0(oct)
write(*,'(a,i5,2g18.10)') 'optcomp,      eps_00(q=0)      :', &
             oct,scrnh0(oct)

        call zinvert_hermitian(scrherm,tm,tmi)
        scri(1:n,1:n)=tmi(:,:)
        ! keep head of inverse dielectric matrix for q=0
        scrnih0(oct)=scri(1,1)

write(*,'(a,i5,2g18.10)') 'optcomp,   eps^-1_00(q=0)      :', &
     oct,scrnih0(oct)
write(*,'(a,i5,2g18.10)') 'optcomp, 1/eps^-1_00(q=0)      :', &
     oct,1.d0/scrnih0(oct)

     end do
     ! symmetrize head of inverse dielectric matrix wrt. the directions
     ! in which q goes to zero
     call symsci0(bsediagsym,scrnh0,scrnih0,scri(1,1))

! set to eps(1,1)
!***scri(1,1)=scrnih0(3)

write(*,'(a,i5,2g18.10)') 'optcomp, symm.   eps^-1_00(q=0):', &
     iqr,scri(1,1)
write(*,'(a,i5,2g18.10)') 'optcomp, symm. 1/eps^-1_00(q=0):', &
          iqr,1.d0/scri(1,1)


!********************************
 ! set upper triangle     
 forall (ig1=2:n,ig2=2:n,ig1.le.ig2) s3(ig1+2,ig2+2)=scrn(ig1,ig2)
!!!        forall (ig1=2:n,ig2=2:n,ig1.gt.ig2) s3(ig1+2,ig2+2)=conjg(s3(ig2+2,ig1+2))
 call zinvert_hermitian(2,s3,s3i)
!!! scri(2:,2:)=s3i(4:,4:)

!****************************************

     ! average the angular dependence of the inverse dielectric matrix at the
     ! Gamma q-point
!!!     call angavdm(n,l,s,bi,scri)
     deallocate(b,bi,u,s)

  else
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
  integer, intent(in) :: n
  complex(8), intent(in) :: bi(n-1,n-1),s(n-1,3),l(3,3)
  complex(8), intent(out) :: av(n,n)
  ! local variables
  integer :: i,j,lm
  real(8) :: t1
  real(8), allocatable :: plat(:,:),p(:)
  complex(8), allocatable :: m00lm(:),m0xlm(:)
  complex(8), allocatable :: ei00(:),ei0x(:,:),ei00lm(:),ei0xlm(:,:)
  allocate(plat(3,lmmaxvr),p(lmmaxvr))
  allocate(m00lm(lmmaxvr),m0xlm(lmmaxvr))
  allocate(ei00(lmmaxvr),ei0x(lmmaxvr,n-1))
  allocate(ei00lm(lmmaxvr),ei0xlm(lmmaxvr,n-1))
  ! unit vectors of spherical covering set in lattice coordinates
  plat=matmul(binv,sphcov)
  ! distances to unit cell boundary in reciprocal space
  p(:)=1.d0/(2.d0*maxval(abs(plat),1))
  do lm=1,lmmaxvr
     ! head which is the 1/(p*L*p) factor
     ei00(lm)=1.d0/dot_product(sphcov(:,lm),matmul(l,sphcov(:,lm)))
     ! wings: -p*s * ei00
     ei0x(lm,:)=-ei00(lm)*matmul(s,sphcov(:,lm))
  end do
! spherical approximation
!***p(:)=((3.d0*(twopi**3/omega))/(fourpi*product(ngridq)))**(1.d0/3.d0)
  ! forward transform shape dependent factor for head to spherical harmonics
  m00lm(:)=matmul(zfshtvr,p)
  ! forward transform shape dependent factor for wings to spherical harmonics
  m0xlm(:)=matmul(zfshtvr,(1.d0/2.d0)*p**2)
  !forward transform to spherical harmonics
  ei00lm=matmul(zfshtvr,ei00)
  ei0xlm=matmul(zfshtvr,ei0x)
  ! angular average
  t1=(omega/(twopi)**3)/product(ngridq)
  av(1,1)=t1*dot_product(m00lm,ei00lm)
  av(1,2:)=t1*matmul(m0xlm,ei0xlm)
  ! body without angular average
  av(2:,2:)=bi(:,:)
  ! set lower triangle
  do i=1,n
     do j=i+1,n
        av(j,i)=conjg(av(i,j))
     end do
  end do 
  deallocate(plat,p,m00lm,m0xlm,ei00,ei0x,ei00lm,ei0xlm)
  
  write(7788,'(2i8,3g18.10)') ((i,j,av(i,j),abs(av(i,j))**2,j=1,n),i=1,n)
   
end subroutine angavdm
