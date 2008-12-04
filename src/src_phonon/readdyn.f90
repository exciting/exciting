
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
integer iq,is,js,ia,ja
integer ip,jp,i,j
real(8) a,b
character(256) fext
! external functions
integer gcd
external gcd
do iq=1,nqpt
  i=0
  do is=1,nspecies
    do ia=1,natoms(is)
      do ip=1,3
        i=i+1
        call phfext(iq,is,ia,ip,fext)
        inquire(file='DYN'//trim(fext),exist=exist)
        if (.not.exist) then
          write(*,*)
          write(*,'("Error(readdyn): file not found :")')
          write(*,'(A)') ' DYN'//trim(fext)
          write(*,*)
          stop
        end if
        open(50,file='DYN'//trim(fext),action='READ',status='OLD', &
         form='FORMATTED')
        j=0
        do js=1,nspecies
          do ja=1,natoms(js)
            do jp=1,3
              j=j+1
              read(50,*) a,b
              dynq(i,j,iq)=cmplx(a,b,8)
            end do
          end do
        end do
        close(50)
      end do
! end loops over atoms and species
    end do
  end do
! symmetrise the dynamical matrix
  call dynsym(vql(:,iq),dynq(:,:,iq))
! end loop over q-vectors
end do
return
end subroutine

