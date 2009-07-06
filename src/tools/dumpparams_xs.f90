
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: dumpparams_xs
! !INTERFACE:
subroutine dumpparams_xs(string,comment)
! !USES:
  use modmain
  use modxs
! !DESCRIPTION:
!   Writes out all input parameters which can be specified in the input file
!   {\tt exciting.in}.
!   Only show those array elements that are within a corresponding cutoff.
!   Trailling whitespaces in string expressions are trimmed.
!   This routine refers only to the parameters related to the excited states
!   implementation.
!
! !REVISION HISTORY:
!   Created July 2008 (Sagmeister)
!EOP
!BOC
  implicit none
  ! arguments
  character(*), intent(in) :: string,comment
  ! local variables
  integer :: i
  call xssetversion
  open(unit=77,file=trim(string),action='write',position='append')
  write(77,*)
  write(77,'("! EXCITING version ",I1.1,".",I1.1,".",I3.3)') version
  write(77,'("! xs (eXited States) version ",I1.1,".",I3.3)') versionxs
  write(77,'(a)') trim(comment)
  write(77,*)
  write(77,'("vgqlmt")')
  write(77,*) nqptmt
  do i=1,nqptmt
    write(77,*) vgqlmt(:,i)
  end do
  write(77,*)
  write(77,'("mdfqtype")')
  write(77,*) mdfqtype
  write(77,*)
  write(77,'("gqmax")')
  write(77,*) gqmax
  write(77,*)
  write(77,'("lmaxapwwf")')
  write(77,*) lmaxapwwf
  write(77,*)
  write(77,'("fastpmat")')
  write(77,*) fastpmat
  write(77,*)
  write(77,'("fastemat")')
  write(77,*) fastemat
  write(77,*)
  write(77,'("emattype")')
  write(77,*) emattype
  write(77,*)
  write(77,'("lmaxemat")')
  write(77,*) lmaxemat
  write(77,*)
  write(77,'("torddf")')
  write(77,*) "'"//trim(torddf)//"'"
  write(77,*)
  write(77,'("tordfxc")')
  write(77,*) "'"//trim(tordfxc)//"'"
  write(77,*)
  write(77,'("acont")')
  write(77,*) acont
  write(77,*)
  write(77,'("nwacont")')
  write(77,*) nwacont
  write(77,*)
  write(77,'("broad")')
  write(77,*) broad
  write(77,*)
  write(77,'("aresdf")')
  write(77,*) aresdf
  write(77,*)
  write(77,'("epsdfde")')
  write(77,*) epsdfde
  write(77,*)
  write(77,'("emaxdf")')
  write(77,*) emaxdf
  write(77,*)
  write(77,'("dfoffdiag")')
  write(77,*) dfoffdiag
  write(77,*)
  write(77,'("tetradf")')
  write(77,*) tetradf
  write(77,*)
  write(77,'("kerndiag")')
  write(77,*) kerndiag
  write(77,*)
  write(77,'("fxctype")')
  write(77,*) fxctype
  write(77,*)
  write(77,'("nexcitmax")')
  write(77,*) nexcitmax
  write(77,*)
  write(77,'("alphalrc")')
  write(77,*) alphalrc
  write(77,*)
  write(77,'("alphalrcdyn")')
  write(77,*) alphalrcdyn
  write(77,*)
  write(77,'("betalrcdyn")')
  write(77,*) betalrcdyn
  write(77,*)
  write(77,'("dftrans")')
  write(77,*) ndftrans
  do i=1,ndftrans
     write(77,*) dftrans(:,i)
  end do
  write(77,*)
  write(77,'("gather")')
  write(77,*) gather
  write(77,*)
  write(77,'("symmorph")')
  write(77,*) symmorph
  write(77,*)
  write(77,'("tevout")')
  write(77,*) tevout
  write(77,*)
  write(77,'("appinfo")')
  write(77,*) tappinfo
  write(77,*)
  write(77,'("dbglev")')
  write(77,*) dbglev
  write(77,*)
  write(77,'("screentype")')
  write(77,*) "'"//trim(screentype)//"'"
  write(77,*)
  write(77,'("nosymscr")')
  write(77,*) nosymscr
  write(77,*)
  write(77,'("reducekscr")')
  write(77,*) reducekscr
  write(77,*)
  write(77,'("ngridkscr")')
  write(77,*) ngridkscr
  write(77,*)
  write(77,'("vkloffscr")')
  write(77,*) vkloffscr
  write(77,*)
  write(77,'("rgkmaxscr")')
  write(77,*) rgkmaxscr
  write(77,*)
  write(77,'("nemptyscr")')
  write(77,*) nemptyscr
  write(77,*)
  write(77,'("scrherm")')
  write(77,*) scrherm
  write(77,*)
  write(77,'("bsetype")')
  write(77,*) "'"//trim(bsetype)//"'"
  write(77,*)
  write(77,'("nosymbse")')
  write(77,*) nosymbse
  write(77,*)
  write(77,'("reducekbse")')
  write(77,*) reducekbse
  write(77,*)
  write(77,'("vkloffbse")')
  write(77,*) vkloffbse
  write(77,*)
  write(77,'("nstlce")')
  write(77,*) nbfce,nafce
  write(77,*)
  write(77,'("nstlbse")')
  write(77,*) nbfbse,nafbse
  close(77)
end subroutine dumpparams_xs
!EOC
