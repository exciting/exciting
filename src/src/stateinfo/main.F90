
! Copyright (C) 2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

program main
  implicit none
  logical :: exis, run
  ! check for binary file
  run=.false.
  inquire(file="STATE.OUT", exist=exis)
  if (exis) then
    call portstate(-1)
    run=.true.
  end if
  ! check for XML file
  inquire(file="STATE.xml", exist=exis)
  if (exis) then
    call portstate(-2)
    run=.true.
  end if
  if (.not. run) then
    write(*,*)
    write(*,'("Error(stateinfo): Neither STATE.OUT nor STATE.xml file found in current working directory.")')
    write(*,*)
  end if
end program

!BOI
! !TITLE: The StateInfo Manual\\ Version 0.7.1
! !AUTHORS: S. Sagmeister and C. Ambrosch-Draxl
! !AFFILIATION:
! !INTRODUCTION: Introduction
! {\sf StateInfo} is a utility which reads in an exciting {\tt STATE.OUT} or
! {\tt STATE.xml} file (the place where the density and the potential as well
! as the magnetization and the B-field are stored) and outputs a short summary
! of its content. In particular, the number of species and atoms, the radial
! meshes, grid sizes, angular momenta cutoffs and spin-polarization are
! reported. It may be used to retrieve some information about the circumstances
! under which the {\tt STATE}-file was produced.
! \section{Usage}
! No input file is required. The program expects either {\tt STATE.OUT} or
! {\tt STATE.xml} (or both) to be located in the current working directory.
!EOI
