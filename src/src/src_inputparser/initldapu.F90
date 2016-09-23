
! Copyright (C) 2011 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine initldapu
      use modinput
      use mod_lda_lu
      Implicit None
      integer :: is
      ldapu=input%groundstate%ldapunumber
      ! zero default values for angular momenta
      llu(:)=0
      ! zero default values for J and U
      ujlu(:,:)=0.d0
      do is=1,size(input%structure%speciesarray)
        if (associated(input%structure%speciesarray(is)%species%LDAplusU)) then
          llu(is)=input%structure%speciesarray(is)%species%LDAplusU%l
          ujlu(1,is)=input%structure%speciesarray(is)%species%LDAplusU%U
          ujlu(2,is)=input%structure%speciesarray(is)%species%LDAplusU%J
        end if
      end do
end subroutine
