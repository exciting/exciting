
! Copyright (C) 2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

program main
  implicit none
  logical :: exis, tbin, txml
  ! check for binary file
  tbin=.false.
  inquire(file="STATE.OUT", exist=exis)
  if (exis) tbin=.true.
  ! check for XML file
  txml=.false.
  inquire(file="STATE.xml", exist=exis)
  if (exis) txml=.true.
  if (tbin .and. txml) then
    write(*,*)
    write(*,'("Error(stateconvert): Both STATE.OUT and STATE.xml files present in current working directory.")')
    write(*,'(" conversion ambiguous")')
    write(*,*)
    stop
  end if
  if ((.not. tbin) .and. (.not. txml)) then
    write(*,*)
    write(*,'("Error(stateconvert): Neither STATE.OUT nor STATE.xml file present in current working directory.")')
    write(*,*)
    stop
  end if
  if (tbin) call portstate(1)
  if (txml) call portstate(2)
end program

!BOI
! !TITLE: The StateConvert Manual\\ Version 0.7.1
! !AUTHORS: S. Sagmeister and C. Ambrosch-Draxl
! !AFFILIATION:
! !INTRODUCTION: Introduction
! {\sf StateConvert} is a utility which converts between the exciting
! {\tt STATE.OUT} binary file and the {\tt STATE.xml} XML file
! (the place where the density and the potential as well
! as the magnetization and the B-field are stored).
! \section{Usage}
! No input file is required. The program expects either {\tt STATE.OUT} or
! {\tt STATE.xml} to be located in the current working directory and performs
! a conversion to the XML format if a binary {\tt STATE}-file is present or to
! the binary format if an XML {\tt STATE}-file is present.
!EOI
