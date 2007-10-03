
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: linengy
! !INTERFACE:
subroutine linengy
! !USES:
use modmain
! !DESCRIPTION:
!   Calculates the new linearisation energies for both the APW and local-orbital
!   radial functions. See the routine {\tt findband}.
!
! !REVISION HISTORY:
!   Created May 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ja,ias,jas
integer l,ilo,io1,io2
! automatic arrays
logical done(natmmax)
real(8) vr(nrmtmax)
! begin loops over atoms and species
do is=1,nspecies
  done(:)=.false.
  do ia=1,natoms(is)
    if (.not.done(ia)) then
      ias=idxas(ia,is)
      vr(1:nrmt(is))=veffmt(1,1:nrmt(is),ias)*y00
!-----------------------!
!     APW functions     !
!-----------------------!
      do l=0,lmaxapw
        do io1=1,apword(l,is)
          if (apwve(io1,l,is)) then
! check if previous radial functions have same default energies
            do io2=1,io1-1
              if (apwve(io2,l,is)) then
                if (abs(apwe0(io1,l,is)-apwe0(io2,l,is)).lt.1.d-4) then
                  apwe(io1,l,ias)=apwe(io2,l,ias)
                  goto 10
                end if
              end if
            end do
! find the band energy starting from default
            apwe(io1,l,ias)=apwe0(io1,l,is)
            call findband(l,0,nprad,nrmt(is),spr(1,is),vr,deband, &
             apwe(io1,l,ias))
          end if
10 continue
        end do
      end do
!---------------------------------!
!     local-orbital functions     !
!---------------------------------!
      do ilo=1,nlorb(is)
        do io1=1,lorbord(ilo,is)
          if (lorbve(io1,ilo,is)) then
! check if previous radial functions have same default energies
            do io2=1,io1-1
              if (lorbve(io2,ilo,is)) then
                if (abs(lorbe0(io1,ilo,is)-lorbe0(io2,ilo,is)).lt.1.d-4) then
                  lorbe(io1,ilo,ias)=lorbe(io2,ilo,ias)
                  goto 20
                end if
              end if
            end do
            l=lorbl(ilo,is)
! find the band energy starting from default
            lorbe(io1,ilo,ias)=lorbe0(io1,ilo,is)
            call findband(l,0,nprad,nrmt(is),spr(1,is),vr,deband, &
             lorbe(io1,ilo,ias))
          end if
20 continue
        end do
      end do
      done(ia)=.true.
! copy to equivalent atoms
      do ja=1,natoms(is)
        if ((.not.done(ja)).and.(eqatoms(ia,ja,is))) then
          jas=idxas(ja,is)
          do l=0,lmaxapw
            do io1=1,apword(l,is)
              apwe(io1,l,jas)=apwe(io1,l,ias)
            end do
          end do
          do ilo=1,nlorb(is)
            do io1=1,lorbord(ilo,is)
              lorbe(io1,ilo,jas)=lorbe(io1,ilo,ias)
            end do
          end do
          done(ja)=.true.
        end if
      end do
    end if
! end loops over atoms and species
  end do
end do
return
end subroutine
!EOC
