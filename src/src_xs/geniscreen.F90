
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine geniscreen(iqr,nmax,n,scri)
  use modxs
  use invert
  implicit none
  ! arguments
  integer, intent(in) :: iqr,n,nmax
  complex(8), intent(out) :: scri(nmax,nmax)
  ! local variables
  complex(8), allocatable :: tm(:,:),tmi(:,:)
  complex(8), allocatable :: scrn(:,:),scrnw(:,:,:),scrnh(:), s3(:,:),s3i(:,:)
  complex(8) :: scrnh0(3),scrnih0(3)
  integer :: oct,oct2,j,ig1,ig2
  integer, external :: octmap

  allocate(scrn(n,n),scrnw(n,2,3),scrnh(9),tm(n,n),tmi(n,n),s3(n+2,n+2),s3i(n+2,n+2))
  ! read screening from file
  call getscreen(iqr,n,scrnh,scrnw,scrn)

  ! invert dielectric matrix for q-point going to zero in x-, y-, and
  ! z-direction for q=0
  if (iqr.eq.1) then
     do oct=1,3
        do oct2=1,3
           j=octmap(oct,oct2)
           s3(oct,oct2)=scrnh(j)
        end do
        if (n.gt.1) then
           s3(oct,4:)=scrnw(2:,1,oct)
        end if
   
        ! index for diagonal tensor component
        j=octmap(oct,oct)
        scrn(1,1)=scrnh(j)
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
