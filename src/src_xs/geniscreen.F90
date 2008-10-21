
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
  implicit none
  ! arguments
  integer, intent(in) :: n
  complex(8), intent(in) :: bi(n-1,n-1),s(n-1,3),l(3,3)
  complex(8), intent(out) :: av(n,n)
  ! local variables
  




end subroutine angavdm
