!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: fermisurf_dx
! !INTERFACE:
!
!
Subroutine fermisurf_dx
! !USES:
      Use modmain
! !DESCRIPTION:
!   Add-on to the {\tt fermisurf} routine.
!   Generates input files for the visualization of the
!   Fermi surface in the first Brillouin zone with OpenDX.
!   Prerequisite is a previous run of {\tt fermisurf}
!   for a $2 \times 2 \times 2$ supercell of reciprocal
!   asymmetric units.
!   It generates a {\tt FERMI.dx} and a {\tt FERMI\_BZ.dx} file
!   where the data for the visualization is stored together
!   with the clipping planes for the construction of the
!   first Brillouin zone.
!   The OpenDX input is complete together with the template
!   {\tt FERMI.net} and the macro {\tt clipvmacro.net}.
!
! !REVISION HISTORY:
!   Created Jan 2005 (Sagmeister)
!   Revised, May 2008 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! local variables
      Integer :: hndl_fermi, hndl_fermidx, hndl_fermidx_bz
      Integer :: nk, nkx, nky, nkz, nst
      Integer :: a, j, ik, obj_startoffset, allocstat, datetime_i (8)
      Integer :: nnlc (26, 3)
      Character (256) :: dump (6)
      Character * 10 file
      Character * 8 date_date
      Character * 10 date_time
      Character * 5 date_zone
      Character * 10, Allocatable :: pos (:), pos_offset (:)
      Real (8) :: nncc (26, 3), b123 (3)
      Real (8), Allocatable :: kvect (:, :), energy (:, :)
!
      Call date_and_time (date_date, date_time, date_zone, datetime_i)
!
  ! initialize
      Call init0
      Call init1
!
  ! the first and second nearest neighbours wrt. to the origin
                       ! end top slice
                       ! end mid slice
                           ! end bottom slice
      nnlc = transpose (reshape( (/ (/ 1, 0, 1 /), (/ 1, 1, 1 /), (/ 0, &
     & 1, 1 /), (/-1, 1, 1 /), (/-1, 0, 1 /), (/-1,-1, 1 /), (/ 0,-1, 1 &
     & /), (/ 1,-1, 1 /), (/ 0, 0, 1 /), (/ 1, 0, 0 /), (/ 1, 1, 0 /), &
     & (/ 0, 1, 0 /), (/-1, 1, 0 /), (/-1, 0, 0 /), (/-1,-1, 0 /), (/ &
     & 0,-1, 0 /), (/ 1,-1, 0 /), (/ 1, 0,-1 /), (/ 1, 1,-1 /), (/ 0, &
     & 1,-1 /), (/-1, 1,-1 /), (/-1, 0,-1 /), (/-1,-1,-1 /), (/ 0,-1,-1 &
     & /), (/ 1,-1,-1 /), (/ 0, 0,-1 /) /), (/ 3, 26 /)))
!
  ! vector to shift the k-grid
      b123 (:) = bvec (:, 1) + bvec (:, 2) + bvec (:, 3)
!
  ! calculate the normal vectors of the clipping planes
      Do a = 1, 26
         nncc (a, :) = nnlc (a, 1) * bvec (:, 1) + nnlc (a, 2) * bvec &
        & (:, 2) + nnlc (a, 3) * bvec (:, 3)
      End Do
      nncc (:, :) = nncc (:, :) * 0.5d0
!
      hndl_fermi = 1101
      hndl_fermidx = 1102
      hndl_fermidx_bz = 1103
      Open (hndl_fermi, File='FERMI.OUT', Action='READ', Form='FORMATTE&
     &D')
      Read (hndl_fermi,*) nkx, nky, nkz, nst, dump
      nk = nkx * nky * nkz
!
      Allocate (kvect(nk, 3), energy(nk, nst), pos(-1:nst), &
     & pos_offset(-1:nst), Stat=allocstat)
      If (allocstat .Ne. 0) Then
         Write (*,*) ' Error(fermisurf_dx): ', 'Memory allocation faile&
        &d.'
         Stop
      End If
!
  ! read in the content of FERMI.OUT
      Do ik = 1, nk
         Read (hndl_fermi,*) (kvect(ik, a), a=1, 3), (energy(ik, a), &
        & a=1, nst)
!!$     write(hndl_fermidx,'(40F12.6)')  (kvect(nk,a),a=1,3), (energy(nk,a), &
!!$          a=1,nst)
      End Do
      Close (hndl_fermi)
!
  ! use this trick with an internal file to convert an integer to a string
      obj_startoffset = 100
      Do a = 1, nst
         Write (File, '(I6)') a - 1
         Read (File, '(A)') pos (a)
         Write (File, '(I6)') obj_startoffset + a + 2
         Read (File, '(A)') pos_offset (a)
      End Do
      Write (File, '(I6)') obj_startoffset + 1
      Read (File, '(A)') pos_offset (-1)
      Write (File, '(I6)') obj_startoffset + 2
      Read (File, '(A)') pos_offset (0)
!
  !
  ! write out to FERMI_BZ.dx
  !
      Write (*,*) ' Writing out to file `FERMI_BZ.dx' // "'" // '...'
      Open (hndl_fermidx_bz, File='FERMI_BZ.dx', Action='WRITE', Form='&
     &FORMATTED')
      Write (hndl_fermidx_bz, '("# Added by EXCITING version ", I2.2, "&
     &.", I2.2, ".", I3.3)') version
      Write (hndl_fermidx_bz, '("# Copyright (C) 2002-2008 J. K. Dewhur&
     &st, S. Sharma and C. Ambrosch-Draxl.")')
      Write (hndl_fermidx_bz, '("# Add-on for plotting Fermisurface wit&
     &h OpenDX)")')
      Write (hndl_fermidx_bz, '("# Copyright (C) 2004-2008 S. Sagmeiste&
     &r and C. Ambrosch-Draxl")')
      Write (hndl_fermidx_bz, '(14A)') '# ', date_date (1:4), ' - ', &
     & date_date (5:6), ' - ', date_date (7:8), '  ', date_time (1:2), &
     & ':', date_time (3:4), ':', date_time (5:6), '  GMT ', date_zone
  ! the normals of the clipping planes for the first Brillouin zone
      Write (hndl_fermidx_bz, '(A)') 'object "clipnormals" class array &
     &type float rank 1 shape 3 items 26 data follows'
      Do a = 1, 26
         Write (hndl_fermidx_bz, '(3F12.6)') (nncc(a, j), j=1, 3)
      End Do
      Write (hndl_fermidx_bz, '(A)') '#'
  ! epilog
      Write (hndl_fermidx_bz, '(A)') ' '
      Write (hndl_fermidx_bz, '(A)') 'end'
      Close (hndl_fermidx_bz)
      Write (*,*) ' ...done.'
!
  !
  ! write out to FERMI.dx
  !
      Write (*,*) ' Writing out to file `FERMI.dx' // "'" // '...'
      Open (hndl_fermidx, File='FERMI.dx', Action='WRITE', Form='FORMAT&
     &TED')
      Write (hndl_fermidx, '("# Added by EXCITING version ", I2.2, ".",&
     & I2.2, ".", I3.3)') version
      Write (hndl_fermidx_bz, '("# Copyright (C) 2002-2008 J. K. Dewhur&
     &st, S. Sharma and C. Ambrosch-Draxl.")')
      Write (hndl_fermidx_bz, '("# Add-on for plotting Fermisurface wit&
     &h OpenDX)")')
      Write (hndl_fermidx_bz, '("# Copyright (C) 2004-2008 S. Sagmeiste&
     &r and C. Ambrosch-Draxl")')
      Write (hndl_fermidx, '(14A)') '# ', date_date (1:4), ' - ', &
     & date_date (5:6), ' - ', date_date (7:8), '  ', date_time (1:2), &
     & ':', date_time (3:4), ':', date_time (5:6), '  GMT ', date_zone
  ! the positions
      Write (hndl_fermidx, '(A, I9, A)') 'object ' // trim &
     & (adjustl(pos_offset(0))) // ' class array type float rank 1 shap&
     &e 3 items ', nk, ' data follows'
      Do ik = 1, nk
         Write (hndl_fermidx, '(3F12.6)') (kvect(ik, a)-b123(a), a=1, &
        & 3)
      End Do
      Write (hndl_fermidx, '(A)') '#'
  ! the connections
      Write (hndl_fermidx, '(A, I5, I5, I5)') 'object ' // trim &
     & (adjustl(pos_offset(-1))) // ' class gridconnections counts ', &
     & nkx, nky, nkz
      Write (hndl_fermidx, '(A)') 'attribute "element type" string "cub&
     &es"'
      Write (hndl_fermidx, '(A)') 'attribute "dep" string "connections"&
     &'
      Write (hndl_fermidx, '(A)') 'attribute "ref" string "positions"'
      Write (hndl_fermidx, '(A)') '#'
  ! the data for the fields
      Do a = 1, nst
         Write (hndl_fermidx, '(A, I9, A)') 'object ' // trim &
        & (adjustl(pos_offset(a))) // ' class array type float rank 0 i&
        &tems ', nk, ' data follows'
         Do ik = 1, nk
            Write (hndl_fermidx, '(F12.6)') energy (ik, a)
         End Do
         Write (hndl_fermidx, '(A)') 'attribute "dep" string "positions&
        &"'
         Write (hndl_fermidx, '(A)') '#'
      End Do
  ! the fields
      Do a = 1, nst
         Write (hndl_fermidx, '(A)') 'object "field' // trim &
        & (adjustl(pos(a))) // '" class field'
         Write (hndl_fermidx, '(A)') 'component "data" value ' // trim &
        & (adjustl(pos_offset(a)))
         Write (hndl_fermidx, '(A)') 'component "positions" value ' // &
        & trim (adjustl(pos_offset(0)))
         Write (hndl_fermidx, '(A)') 'component "connections" value ' &
        & // trim (adjustl(pos_offset(-1)))
         Write (hndl_fermidx, '(A)') 'attribute "name" string "field' &
        & // trim (adjustl(pos(a))) // '"'
         Write (hndl_fermidx, '(A)') '#'
      End Do
  ! the default object
      Write (hndl_fermidx, '(A)') 'object "default" class group'
      Do a = 1, nst
         Write (hndl_fermidx, '(A)') 'member "field' // trim &
        & (adjustl(pos(a))) // '" value "field' // trim &
        & (adjustl(pos(a))) // '"'
      End Do
      Write (hndl_fermidx, '(A)') '#'
  ! epilog
      Write (hndl_fermidx, '(A)') ' '
      Write (hndl_fermidx, '(A)') 'end'
      Close (hndl_fermidx)
      Write (*,*) ' ...done.'
!
  ! deallocate
      Deallocate (kvect, energy, pos, pos_offset)
!
End Subroutine fermisurf_dx
!EOC
