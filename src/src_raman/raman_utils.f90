! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created July 2014 by Stefan Kontur
!
!
Module m_raman_utils
      Implicit None
Contains
!
!
   Subroutine raman_save_struct
      Use mod_atoms
      use mod_lattice, only: ainv
      Use modinput
      use mod_phonon, only: natoms0, natmtot0, avec0, ainv0, atposc0
      Implicit None
      integer :: is, ia
      natoms0 (1:nspecies) = natoms (1:nspecies)
      natmtot0 = natmtot
      avec0 (:, :) = input%structure%crystal%basevect(:, :)
      ainv0 (:, :) = ainv (:, :)
      atposc0 (:, :, :) = 0.d0
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            atposc0 (:, ia, is) = atposc (:, ia, is)
         End Do
      End Do
   end subroutine raman_save_struct
!
!
   Subroutine raman_restore_struct
      Use mod_atoms
      use mod_lattice, only: ainv
      Use modinput
      use mod_phonon, only: natoms0, natmtot0, avec0, ainv0, atposc0
      Implicit None
      integer :: is, ia
      natoms (1:nspecies) = natoms0(1:nspecies)
      natmtot = natmtot0
      input%structure%crystal%basevect(:, :) = avec0(:, :)
      ainv (:, :) = ainv0 (:, :)
      atposc (:, :, :) = 0.d0
      Do is = 1, nspecies
       Do ia = 1, natoms (is)
        ! restore cartesian coordinates
        atposc (:, ia, is) = atposc0 (:, ia, is)
        ! restore lattice coordinates in input structure
        call r3mv (ainv, atposc(:, ia, is), input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:))
       End Do
      End Do
   end subroutine raman_restore_struct
!
!
   subroutine raman_readpot (istep, fnam, dph, engy, force)
      use m_getunit
      implicit none
    ! arguments
      Character (*), Intent (In) :: fnam
      integer, intent(out) :: istep
      real(8), intent(out) :: dph, engy, force
    ! local variables
      logical :: existent, opened
      integer :: funit
    ! check if file exists
      Inquire (File=trim(fnam), Exist=existent)
      if (.not. existent) then
        write(*, '("Error(Raman_readpot): File ",a," does not exist!")') trim(fnam)
        stop
      endif
    ! check if file is opened
      Inquire (File=trim(fnam), Opened=Opened, Number=funit)
      If (opened) Then
         Close (funit)
      End If
      Call getunit (funit)
      Open (funit, File=trim(fnam), Action='read')
      read(funit, *) istep, dph, engy, force
      close (funit)
   end subroutine raman_readpot
!
!
   subroutine raman_writepot (istep, fnam, dph, engy, force)
      use m_getunit
      implicit none
    ! arguments
      Character (*), Intent (In) :: fnam
      integer, intent(in) :: istep
      real(8), intent(in) :: dph, engy, force
    ! local variables
      integer :: funit
    !
      call getunit (funit)
      open (unit=funit, file=trim(fnam), status='unknown', action='write')
      write (funit, '(i6,f14.6,g23.15,g21.13)') istep, dph, engy, force
      close (unit=funit)
   end subroutine raman_writepot
!
!
   Subroutine getfgew ( eigv )
    ! computes the factor f relating the atomic displacements u to the normal coordinate Q
    ! it is also used to normalize u and represents the effective mass of the phonon mode
      use raman_ew, only: fgew
      use mod_atoms, only: natmtot, nspecies, natoms, spmass
      implicit none
    ! arguments
      complex(8), intent(in) :: eigv(3*natmtot)  
    ! local variables
      integer :: i, ia, is, iat
    !
      fgew = 0.d0
      iat = 0
      do is = 1, nspecies
         do ia = 1, natoms(is)
            iat = iat + 1
            do i = 1, 3
               fgew = fgew + dble(eigv(3*(iat-1)+i))**2/spmass(is)
            enddo
         enddo
      enddo
      fgew = 1.d0 / fgew
   End Subroutine getfgew
!
!
   subroutine raman_readeps (imode, istep, i1, i2, fnam)
   ! read dielectric function from file
      use modxs
      use raman_coeff
      use m_getunit
      use modinput
      implicit none
    ! arguments
      integer, intent(in) :: imode, istep, i1, i2
      Character (*), Intent (In) :: fnam
    ! local variables
      logical :: existent, opened
      integer :: iw, funit, ncommentlines
      real (8) :: w, eps_re, eps_im, eps_rekk
    ! Set number of comment lines to skip in EPSILON*.OUT
      ncommentlines=18
    ! check if file exists
      Inquire (File=trim(fnam), Exist=existent)
      if (.not. existent) then
        write(*,*) 'File ',trim(fnam),' does not exist!'
        stop
      endif
    ! check if file is opened
      Inquire (File=trim(fnam), Opened=Opened, Number=funit)
      If (opened) Then
         Close (funit)
      End If
      Call getunit (funit)
      Open (funit, File=trim(fnam), Action='read')
      do iw = 1, ncommentlines
         read(funit,*)
      enddo
      do iw = 1, nwdf
         read(funit, *) w, eps_re, eps_im
         df(imode, istep, i1, i2, iw) = cmplx (eps_re, eps_im, 8)
      enddo
      close (funit)
   end subroutine raman_readeps 
!
!
   subroutine raman_writeeps (imode, istep, oct1, oct2, comp_ch, fxt)
      use modxs
      use raman_coeff
      use modinput
      use m_getunit
      implicit none
    ! arguments
      integer, intent(in) :: imode, istep, oct1, oct2
      character(2), intent(in) :: comp_ch
      character(*), intent(in) :: fxt
      real (8) :: o_mega, t1
      integer :: funit, iw
 !
      call getunit (funit)
      open (unit=funit, file='../EPSILON_OC'//comp_ch//trim(fxt), status='unknown', action='write')
      t1 = (input%xs%energywindow%intv(2) - input%xs%energywindow%intv(1)) / dble(input%xs%energywindow%points)
      do iw = 1, nwdf
         o_mega = t1 * (iw - 1) + input%xs%energywindow%intv(1)
         write (funit, '(3g18.10)') o_mega*escale, df(imode, istep, oct1, oct2, iw)
      enddo
      close (unit=funit)
   end subroutine raman_writeeps
!
!
End Module m_raman_utils
