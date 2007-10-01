
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readdyn(dynq)
use modmain
implicit none
! arguments
complex(8), intent(out) :: dynq(3*natmtot,3*natmtot,nqpt)
! local variables
logical exist
integer iq,is,js,ia,ja,ip,jp
integer i,j,m(3),n(3)
real(8) a,b
! external functions
integer gcd
external gcd
do iq=1,nqpt
  do i=1,3
    if (ivq(i,iq).ne.0) then
      j=gcd(ivq(i,iq),ngridq(i))
      m(i)=ivq(i,iq)/j
      n(i)=ngridq(i)/j
    else
      m(i)=0
      n(i)=0
    end if
  end do
  i=0
  do is=1,nspecies
    do ia=1,natoms(is)
      write(filext,'("_Q",2I2.2,"_",2I2.2,"_",2I2.2,"_S",I2.2,"_A",I3.3,&
       &".OUT")') m(1),n(1),m(2),n(2),m(3),n(3),is,ia
      inquire(file='DYN'//trim(filext),exist=exist)
      if (.not.exist) then
        write(*,*)
        write(*,'("Error(readdyn): file not found :")')
        write(*,'(A)') ' DYN'//trim(filext)
        write(*,*)
        stop
      end if
      open(50,file='DYN'//trim(filext),action='READ',status='OLD', &
       form='FORMATTED')
      do ip=1,3
        i=i+1
        read(50,*)
        j=0
        do js=1,nspecies
          do ja=1,natoms(js)
            do jp=1,3
              read(50,*) a,b
              dynq(i,j+jp,iq)=cmplx(a,b,8)
            end do
            j=j+3
          end do
        end do
      end do
      close(50)
! end loops over atoms and species
    end do
  end do
! symmetrise the dynamical matrix
  call dynsym(vql(1,iq),dynq(1,1,iq))
! end loop over q-vectors
end do
return
end subroutine

