!
Module m_raman_utils
      Implicit None
Contains
!
!
   Subroutine raman_delgndst
      Use modmain
      Implicit None
! delete the eigenvector files
      Call delevec
! delete the eigenvalue files
      Open (70, File=trim(scrpath)//'EVALFV'//trim(filext))
      Close (70, Status='DELETE')
      Open (70, File=trim(scrpath)//'EVALSV'//trim(filext))
      Close (70, Status='DELETE')
! delete the occupancy file
      Open (70, File=trim(scrpath)//'OCCSV'//trim(filext))
      Close (70, Status='DELETE')
! delete the STATE.OUT file
      Open (50, File='STATE'//trim(filext))
      Close (50, Status='DELETE')
      Return
   End Subroutine raman_delgndst
!
   subroutine raman_readpot (istep, fnam, dph, engy)
      use m_getunit
      implicit none
    ! arguments
      Character (*), Intent (In) :: fnam
      integer, intent(out) :: istep
      real(8), intent(out) :: dph, engy
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
      read(funit, *) istep, dph, engy
      close (funit)
   end subroutine raman_readpot
!
   subroutine raman_writepot (istep, fnam, dph, engy)
      use m_getunit
      implicit none
    ! arguments
      Character (*), Intent (In) :: fnam
      integer, intent(in) :: istep
      real(8), intent(in) :: dph, engy
    ! local variables
      integer :: funit
    !
      call getunit (funit)
      open (unit=funit, file=trim(fnam), status='unknown', action='write')
      write (funit, '(i6,g18.8,g28.18)') istep, dph, engy
      close (unit=funit)
   end subroutine raman_writepot
!
   subroutine raman_readeps (imode, istep, i1, i2, fnam)
   ! read dielectric function from file
      use modxs
      use raman_coeff
      use raman_input
      use m_getunit
      use modinput
      implicit none
    ! arguments
      integer, intent(in) :: imode, istep, i1, i2
      Character (*), Intent (In) :: fnam
    ! local variables
      logical :: existent, opened
      integer :: iw, funit
      real (8) :: w, eps_re, eps_im, eps_rekk
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
   Subroutine raman_delxs (fnam)
      Use m_getunit
      Implicit None
    ! arguments
      Character (*), Intent (In) :: fnam
    ! local variables
      Integer, Parameter :: verb = 0
      Integer :: funit
      Logical :: existent, opened
    ! check if file exists
      Inquire (File=trim(fnam), Exist=existent)
      If ((verb .Gt. 0) .And. ( .Not. existent)) Then
         Write (*, '("Warning(Raman_delxs): attempted to delete non-existent file: ", a)') trim (fnam)
         Return
      End If
    ! check if file is opened
      Inquire (File=trim(fnam), Opened=Opened, Number=funit)
    ! close file if opened
      If (opened) Then
         Close (funit)
      End If
    ! open file for writing
      Call getunit (funit)
      Open (funit, File=trim(fnam), Action='write')
    ! delete file
      Close (funit, Status='delete')
   End Subroutine raman_delxs
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
      write(*,*) 'fgew ',fgew
   End Subroutine getfgew
!
End Module m_raman_utils
