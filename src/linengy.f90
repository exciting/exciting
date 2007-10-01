
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
integer is,ia1,ia2,ias1,ias2
integer l,ilo,io1,io2
! automatic arrays
logical done(natmmax)
real(8) vr(nrmtmax)
! begin loops over atoms and species
do is=1,nspecies
  done(:)=.false.
  do ia1=1,natoms(is)
    if (.not.done(ia1)) then
      ias1=idxas(ia1,is)
      vr(1:nrmt(is))=veffmt(1,1:nrmt(is),ias1)*y00
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
                  apwe(io1,l,ias1)=apwe(io2,l,ias1)
                  goto 10
                end if
              end if
            end do
! find the band energy starting from default
            apwe(io1,l,ias1)=apwe0(io1,l,is)
            call findband(l,nprad,nrmt(is),spr(1,is),vr,deband,apwe(io1,l,ias1))
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
                  lorbe(io1,ilo,ias1)=lorbe(io2,ilo,ias1)
                  goto 20
                end if
              end if
            end do
            l=lorbl(ilo,is)
! find the band energy starting from default
            lorbe(io1,ilo,ias1)=lorbe0(io1,ilo,is)
            call findband(l,nprad,nrmt(is),spr(1,is),vr,deband, &
             lorbe(io1,ilo,ias1))
          end if
20 continue
        end do
      end do
      done(ia1)=.true.
! copy to equivalent atoms
      do ia2=1,natoms(is)
        if ((.not.done(ia2)).and.(nsymeqat(ia2,ia1,is).gt.0)) then
          ias2=idxas(ia2,is)
          do l=0,lmaxapw
            do io1=1,apword(l,is)
              apwe(io1,l,ias2)=apwe(io1,l,ias1)
            end do
          end do
          do ilo=1,nlorb(is)
            do io1=1,lorbord(ilo,is)
              lorbe(io1,ilo,ias2)=lorbe(io1,ilo,ias1)
            end do
          end do
          done(ia2)=.true.
        end if
      end do
    end if
! end loops over atoms and species
  end do
end do
return
end subroutine
!EOC
