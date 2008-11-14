
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: dumpparams
! !INTERFACE:
subroutine dumpparams(string,comment,sppath_,sc_,sc1_,sc2_,sc3_,vacuum_)
! !USES:
  use modmain
! !DESCRIPTION:
!   Writes out all input parameters which can be specified in the input file
!   {\tt exciting.in}.
!   Only show those array elements that are within a corresponding cutoff.
!   Trailling whitespaces in string expressions are trimmed.
!   This routine refers only to parameters of the main version of Exciting.
!
! !REVISION HISTORY:
!   Created 2007 (Sagmeister)
!   Added parameters, July 2008 (Sagmeister)
!EOP
!BOC
  implicit none
  ! arguments
  character(*), intent(in) :: string,comment
  character(*), intent(in) :: sppath_
  real(8), intent(in) :: sc_,sc1_,sc2_,sc3_
  real(8), intent(in) :: vacuum_
  ! local variables
  integer :: j,ia,is
  open(unit=77,file=trim(string),action='write',status='replace')
  write(77,*)
  write(77,'("! EXCITING version ",I1.1,".",I1.1,".",I3.3)') version
  write(77,'(a)') trim(comment)
  write(77,*)
  write(77,'("tasks")')
  do j=1,ntasks
     write(77,*) tasks(j)
  end do
  write(77,*)
  write(77,'("avec")')
  write(77,*) avec(1,1),avec(2,1),avec(3,1)
  write(77,*) avec(1,2),avec(2,2),avec(3,2)
  write(77,*) avec(1,3),avec(2,3),avec(3,3)
  write(77,*)
  write(77,'("scale")')
  write(77,*) sc_
  write(77,*)
  write(77,'("scale1")')
  write(77,*) sc1_
  write(77,*)
  write(77,'("scale2")')
  write(77,*) sc2_
  write(77,*)
  write(77,'("scale3")')
  write(77,*) sc3_
  write(77,*)
  write(77,'("epslat")')
  write(77,*) epslat
  write(77,*)
  write(77,'("primcell")')
  write(77,*) primcell
  write(77,*)
  write(77,'("tshift")')
  write(77,*) tshift
  write(77,*)
  write(77,'("autokpt")')
  write(77,*) autokpt
  write(77,*)
  write(77,'("ngridk")')
  write(77,*) ngridk(1),ngridk(2),ngridk(3)
  write(77,*)
  write(77,'("vkloff")')
  write(77,*) vkloff(1),vkloff(2),vkloff(3)
  write(77,*)
  write(77,'("reducek")')
  write(77,*) reducek
  write(77,*)
  write(77,'("ngridq")')
  write(77,*) ngridq(1),ngridq(2),ngridq(3)
  write(77,*)
  write(77,'("reduceq")')
  write(77,*) reduceq
  write(77,*)
  write(77,'("rgkmax")')
  write(77,*) rgkmax
  write(77,*)
  write(77,'("gmaxvr")')
  write(77,*) gmaxvr
  write(77,*)
  write(77,'("lmaxapw")')
  write(77,*) lmaxapw
  write(77,*)
  write(77,'("lmaxvr")')
  write(77,*) lmaxvr
  write(77,*)
  write(77,'("lmaxmat")')
  write(77,*) lmaxmat
  write(77,*)
  write(77,'("lmaxinr")')
  write(77,*) lmaxinr
  write(77,*)
  write(77,'("fracinr")')
  write(77,*) fracinr
  write(77,*)
  write(77,'("npsden")')
  write(77,*) npsden
  write(77,*)
  write(77,'("spinpol")')
  write(77,*) spinpol
  write(77,*)
  write(77,'("spinorb")')
  write(77,*) spinorb
  write(77,*)
  write(77,'("xctype")')
  write(77,*) xctype
  write(77,*)
  write(77,'("stype")')
  write(77,*) stype
  write(77,*)
  write(77,'("swidth")')
  write(77,*) swidth
  write(77,*)
  write(77,'("epsocc")')
  write(77,*) epsocc
  write(77,*)
  write(77,'("epschg")')
  write(77,*) epschg
  write(77,*)
  write(77,'("nempty")')
  write(77,*) nempty
  write(77,*)
  write(77,'("beta0")')
  write(77,*) beta0
  write(77,*)
  write(77,'("maxscl")')
  write(77,*) maxscl
  write(77,*)
  write(77,'("epspot")')
  write(77,*) epspot
  write(77,*)
  write(77,'("epsengy")')
  write(77,*) epsengy
  write(77,*)
  write(77,'("epsforce")')
  write(77,*) epsforce
  write(77,*)
  write(77,'("cfdamp")')
  write(77,*) cfdamp
  write(77,*)
  write(77,'("sppath")')
  write(77,*) "'"//trim(sppath_)//"'"
  write(77,*)
  write(77,'("scrpath")')
  write(77,*) "'"//trim(scrpath)//"'"
  write(77,*)
  write(77,'("molecule")')
  write(77,*) molecule
  write(77,*)
  write(77,'("vacuum")')
  write(77,*) vacuum_
  write(77,*)
  write(77,'("atoms")')
  write(77,*) nspecies
  do is=1,nspecies
     write(77,*) "'"//trim(spfname(is))//"'"
     write(77,*) natoms(is)
     do ia=1,natoms(is)
        write(77,*) atposl(1,ia,is),atposl(2,ia,is),atposl(3,ia,is), &
             bfcmt(1,ia,is),bfcmt(2,ia,is),bfcmt(3,ia,is)
     end do
  end do
  write(77,*)
  write(77,'("plot1d")')
  write(77,*) nvp1d,npp1d
  do j=1,nvp1d
     write(77,*) vvlp1d(1,j),vvlp1d(2,j),vvlp1d(3,j)
  end do
  write(77,*)
  write(77,'("plot2d")')
  write(77,*) vclp2d(1,1),vclp2d(2,1),vclp2d(3,1)
  write(77,*) vclp2d(1,2),vclp2d(2,2),vclp2d(3,2)
  write(77,*) vclp2d(1,3),vclp2d(2,3),vclp2d(3,3)
  write(77,*) np2d(1),np2d(2)
  write(77,*)

  write(77,'("dos")')
  write(77,*) nwdos,ngrdos,nsmdos
  write(77,*) wdos(1),wdos(2)
  write(77,*)
  write(77,'("tau0atm")')
  write(77,*) tau0atm
  write(77,*)
  write(77,'("nstfsp")')
  write(77,*) nstfsp
  write(77,*)
  write(77,'("lradstp")')
  write(77,*) lradstp
  write(77,*)
  write(77,'("chgexs")')
  write(77,*) chgexs
  write(77,*)
  write(77,'("nprad")')
  write(77,*) nprad
  write(77,*)
  write(77,'("scissor")')
  write(77,*) scissor
  write(77,*)
  write(77,'("optcomp")')
  do j=1,noptcomp
     write(77,*) optcomp(:,j)
  end do
  write(77,*)
  write(77,'("usegdft")')
  write(77,*) usegdft
  write(77,*)
  write(77,'("intraband")')
  write(77,*) intraband
  write(77,*)
  write(77,'("evalmin")')
  write(77,*) evalmin
  write(77,*)
  write(77,'("deband")')
  write(77,*) deband
  write(77,*)
  write(77,'("bfieldc")')
  write(77,*) bfieldc
  write(77,*)
  write(77,'("fixspin")')
  write(77,*) fixspin
  write(77,*)
  write(77,'("momfix")')
  write(77,*) momfix
  write(77,*)
  write(77,'("mommtfix")')
  do is=1,nspecies
     do ia=1,natoms(is)
        write(77,*) is,ia,mommtfix(:,ia,is)
     end do
  end do
  write(77,*)
  write(77,'("taufsm")')
  write(77,*) taufsm
  write(77,*)
  write(77,'("autormt")')
  write(77,*) autormt
  write(77,*)
  write(77,'("rmtapm")')
  write(77,*) rmtapm
  write(77,*)
  write(77,'("nosym")')
  write(77,*) nosym
  write(77,*)
  write(77,'("deltaph")')
  write(77,*) deltaph
  write(77,*)
  write(77,'("phwrite")')
  write(77,*) nphwrt
  do j=1,nphwrt
     write(77,*) vqlwrt(:,j)
  end do
  write(77,*)
  write(77,'("notes")')
  do j=1,notelns
     write(77,'(A80)') notes(j)
  end do
  write(77,*)
  write(77,'("tforce")')
  write(77,*) tforce
  write(77,*)
  write(77,'("tfibs")')
  write(77,*) tfibs
  write(77,*)
  write(77,'("maxitoep")')
  write(77,*) maxitoep
  write(77,*)
  write(77,'("tauoep")')
  write(77,*) tauoep
  write(77,*)
  write(77,'("kstlist")')
  do j=1,nkstlist
     write(77,*) kstlist(:,j)
  end do
  write(77,*)
  write(77,'("vklem")')
  write(77,*) vklem
  write(77,*)
  write(77,'("deltaem")')
  write(77,*) deltaem
  write(77,*)
  write(77,'("ndspem")')
  write(77,*) ndspem
  write(77,*)
  write(77,'("nosource")')
  write(77,*) nosource
  write(77,*)
  write(77,'("spinsprl")')
  write(77,*) spinsprl
  write(77,*)
  write(77,'("vqlss")')
  write(77,*) vqlss
  write(77,*)
  write(77,'("nwrite")')
  write(77,*) nwrite
  write(77,*)
  write(77,'("tevecsv")')
  write(77,*) tevecsv
  write(77,*)
  write(77,'("lda+u")')
  write(77,*) ldapu
  do is=1,nspecies
     write(77,*) is,llu(is),ujlu(1,is),ujlu(2,is)
  end do
  write(77,*)
  write(77,'("rdmxctype")')
  write(77,*) rdmxctype
  write(77,*)
  write(77,'("rdmmaxscl")')
  write(77,*) rdmmaxscl
  write(77,*)
  write(77,'("maxitn")')
  write(77,*) maxitn
  write(77,*)
  write(77,'("maxitc")')
  write(77,*) maxitc
  write(77,*)
  write(77,'("taurdmn")')
  write(77,*) taurdmn
  write(77,*)
  write(77,'("taurdmc")')
  write(77,*) taurdmc
  write(77,*)
  write(77,'("rdmalpha")')
  write(77,*) rdmalpha
  write(77,*)
  write(77,'("reducebf")')
  write(77,*) reducebf
  close(77)
end subroutine dumpparams
!EOC
