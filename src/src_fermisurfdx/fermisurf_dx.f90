
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: fermisurf_dx
! !INTERFACE:
subroutine fermisurf_dx
! !USES:
use modmain
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
  implicit none
  ! local variables
  integer :: hndl_fermi, hndl_fermidx, hndl_fermidx_bz
  integer :: nk, nkx, nky, nkz, nst
  integer :: a, j, ik, ist, obj_startoffset, allocstat, datetime_i(8)
  integer :: nnlc(26,3)
  character(256) :: dump(6)
  character*10 file
  character*8 date_date
  character*10 date_time
  character*5 date_zone
  character*10, allocatable :: pos(:), pos_offset(:)
  real(8) :: nncc(26,3), b123(3)
  real(8), allocatable :: kvect(:,:), energy(:,:)

  call date_and_time(date_date, date_time, date_zone, datetime_i)

  ! initialize
  call init0
  call init1

  ! the first and second nearest neighbours wrt. to the origin
  nnlc = transpose(reshape( (/ &
       (/1, 0, 1/), &
       (/1, 1, 1/), &
       (/0, 1, 1/), &
       (/-1, 1, 1/), &
       (/-1, 0, 1/), &
       (/-1, -1, 1/), &
       (/0, -1, 1/), &
       (/1, -1, 1/), &
       (/0, 0, 1/), &  ! end top slice
       (/1, 0, 0/), &
       (/1, 1, 0/), &
       (/0, 1, 0/), &
       (/-1, 1, 0/), &
       (/-1, 0, 0/), &
       (/-1, -1, 0/), &
       (/0, -1, 0/), &
       (/1, -1, 0/), &  ! end mid slice
       (/1, 0, -1/), &
       (/1, 1, -1/), &
       (/0, 1, -1/), &
       (/-1, 1, -1/), &
       (/-1, 0, -1/), &
       (/-1, -1, -1/), &
       (/0, -1, -1/), &
       (/1, -1, -1/), &
       (/0, 0, -1/) /),  & ! end bottom slice
       (/3, 26/) ) )

  ! vector to shift the k-grid
  b123(:) = bvec(:,1) + bvec(:,2) + bvec(:,3)

  ! calculate the normal vectors of the clipping planes
  do a = 1, 26
     nncc(a,:) = nnlc(a,1)*bvec(:,1) + nnlc(a,2)*bvec(:,2) + nnlc(a,3)*bvec(:,3)
  end do
  nncc(:,:) = nncc(:,:)*0.5d0

  hndl_fermi = 1101
  hndl_fermidx = 1102
  hndl_fermidx_bz = 1103
  open(hndl_fermi,file='FERMI.OUT',action='READ',form='FORMATTED')
  read(hndl_fermi,*) nkx, nky, nkz, nst, dump
  nk = nkx * nky * nkz

  allocate( kvect(nk,3), &
       energy(nk,nst), &
       pos(-1:nst), &
       pos_offset(-1:nst), &
       STAT = allocstat )
  if ( allocstat .ne. 0 ) then
     write(*,*) ' Error(fermisurf_dx): ', 'Memory allocation failed.'
     stop
  end if

  ! read in the content of FERMI.OUT
  do ik = 1, nk
     read(hndl_fermi, *) (kvect(ik,a),a=1,3), (energy(ik,a),a=1,nst)
!!$     write(hndl_fermidx,'(40F12.6)')  (kvect(nk,a),a=1,3), (energy(nk,a), &
!!$          a=1,nst)
  end do
  close(hndl_fermi)

  ! use this trick with an internal file to convert an integer to a string
  obj_startoffset = 100
  do a = 1, nst
     write(file, '(I6)') a - 1
     read(file, '(A)') pos(a)
     write(file, '(I6)') obj_startoffset + a + 2
     read(file, '(A)') pos_offset(a)
  end do
  write(file, '(I6)') obj_startoffset + 1
  read(file, '(A)') pos_offset(-1)
  write(file, '(I6)') obj_startoffset + 2
  read(file, '(A)') pos_offset(0)

  !
  ! write out to FERMI_BZ.dx
  !
  write(*,*) ' Writing out to file `FERMI_BZ.dx'//"'"//'...'
  open(hndl_fermidx_bz,file='FERMI_BZ.dx',action='WRITE',form='FORMATTED')
  write(hndl_fermidx_bz,'("# Added by EXCITING version ",I2.2,".",I2.2,".",&
       &I3.3)') version
  write(hndl_fermidx_bz,'("# Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma &
       &and C. Ambrosch-Draxl.")')
  write(hndl_fermidx_bz,'("# Add-on for plotting Fermisurface with OpenDX)")')
  write(hndl_fermidx_bz,'("# Copyright (C) 2004-2008 S. Sagmeister and C. &
       &Ambrosch-Draxl")')
  write(hndl_fermidx_bz,'(14A)') '# ',date_date(1:4),'-', &
       date_date(5:6),'-',date_date(7:8),'  ', &
       date_time(1:2),':',date_time(3:4),':',date_time(5:6), &
       '  GMT ', date_zone
  ! the normals of the clipping planes for the first Brillouin zone
  write(hndl_fermidx_bz,'(A)') 'object "clipnormals" class array type float &
       &rank 1 shape 3 items 26 data follows'
  do a = 1, 26
     write(hndl_fermidx_bz,'(3F12.6)') (nncc(a,j),j=1,3)
  end do
  write(hndl_fermidx_bz,'(A)') '#'
  ! epilog
  write(hndl_fermidx_bz,'(A)') ' '
  write(hndl_fermidx_bz,'(A)') 'end'
  close(hndl_fermidx_bz)
  write(*,*) ' ...done.'

  !
  ! write out to FERMI.dx
  !
  write(*,*) ' Writing out to file `FERMI.dx'//"'"//'...'
  open(hndl_fermidx,file='FERMI.dx',action='WRITE',form='FORMATTED')
  write(hndl_fermidx,'("# Added by EXCITING version ",I2.2,".",I2.2,".",&
       &I3.3)') version
  write(hndl_fermidx_bz,'("# Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma &
       &and C. Ambrosch-Draxl.")')
  write(hndl_fermidx_bz,'("# Add-on for plotting Fermisurface with OpenDX)")')
  write(hndl_fermidx_bz,'("# Copyright (C) 2004-2008 S. Sagmeister and C. &
       &Ambrosch-Draxl")')
  write(hndl_fermidx,'(14A)') '# ',date_date(1:4),'-', &
       date_date(5:6),'-',date_date(7:8),'  ', &
       date_time(1:2),':',date_time(3:4),':',date_time(5:6), &
       '  GMT ', date_zone
  ! the positions
  write(hndl_fermidx,'(A,I9,A)') 'object '//trim(adjustl(pos_offset(0)))// &
       ' class array type float rank 1 shape 3 items ', nk, ' data follows'
  do ik = 1, nk
     write(hndl_fermidx,'(3F12.6)') (kvect(ik,a)-b123(a),a=1,3)
  end do
  write(hndl_fermidx,'(A)') '#'
  ! the connections
  write(hndl_fermidx,'(A,I5,I5,I5)') 'object '//trim(adjustl(pos_offset(-1)))//&
       ' class gridconnections counts ', nkx, nky, nkz
  write(hndl_fermidx,'(A)') 'attribute "element type" string "cubes"'
  write(hndl_fermidx,'(A)') 'attribute "dep" string "connections"'
  write(hndl_fermidx,'(A)') 'attribute "ref" string "positions"'
  write(hndl_fermidx,'(A)') '#'
  ! the data for the fields
  do a = 1, nst
     write(hndl_fermidx,'(A,I9,A)') 'object '//trim(adjustl(pos_offset(a)))// &
          ' class array type float rank 0 items ', nk, ' data follows'
     do ik = 1, nk
        write(hndl_fermidx,'(F12.6)') energy(ik,a)
     end do
     write(hndl_fermidx,'(A)') 'attribute "dep" string "positions"' 
     write(hndl_fermidx,'(A)') '#'
  end do
  ! the fields
  do a = 1, nst
     write(hndl_fermidx,'(A)') 'object "field'//trim(adjustl(pos(a)))// &
          '" class field'
     write(hndl_fermidx,'(A)') 'component "data" value '// &
          trim(adjustl(pos_offset(a)))
     write(hndl_fermidx,'(A)') 'component "positions" value '// &
          trim(adjustl(pos_offset(0)))
     write(hndl_fermidx,'(A)') 'component "connections" value '// &
          trim(adjustl(pos_offset(-1)))
     write(hndl_fermidx,'(A)') 'attribute "name" string "field'// &
          trim(adjustl(pos(a)))//'"'
     write(hndl_fermidx,'(A)') '#'
  end do
  ! the default object
  write(hndl_fermidx,'(A)') 'object "default" class group'
  do a = 1, nst
     write(hndl_fermidx,'(A)') 'member "field'//trim(adjustl(pos(a)))// &
          '" value "field'//trim(adjustl(pos(a)))//'"'
  end do
  write(hndl_fermidx,'(A)') '#'
  ! epilog
  write(hndl_fermidx,'(A)') ' '
  write(hndl_fermidx,'(A)') 'end'
  close(hndl_fermidx)
  write(*,*) ' ...done.'

  ! deallocate
  deallocate(kvect, energy, pos, pos_offset)

end subroutine fermisurf_dx
!EOC
