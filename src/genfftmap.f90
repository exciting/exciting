!
!BOP
! !ROUTINE: genfftmap
! !INTERFACE:
!
!
Subroutine genfftmap(fftmap,gmaxcustom)
! !USES:
      Use modinput
      Use modmain
      use modxs, only : fftmap_type
! !DESCRIPTION:
! Prepares a new FFT grid with given $G_\mathrm{max}$.
!  
! !REVISION HISTORY:
!   Created October 2014 (Andris)
!EOP
!BOC
      Implicit None
      type (fftmap_type) :: fftmap
      integer :: intgv_custom (3, 2)
      integer :: i1,i2,i3,ig,j1,j2,j3
      real(8) :: gmaxcustom
      
! find optimal grid size for potential and density
      fftmap%ngrid (:) = Int (gmaxcustom*&
     & Sqrt(input%structure%crystal%basevect(1, :)**2+&
     & input%structure%crystal%basevect(2, :)**2+&
     & input%structure%crystal%basevect(3, :)**2)/pi) + 1
#ifndef FFTW
! find next largest FFT-compatible grid size
      Call nfftifc (fftmap%ngrid(1))
      Call nfftifc (fftmap%ngrid(2))
      Call nfftifc (fftmap%ngrid(3))
#endif

      If ((fftmap%ngrid(1) .Le. 0) .Or. (fftmap%ngrid(2) .Le. 0) .Or. (fftmap%ngrid(3) .Le. &
     & 0)) Then
         Write (*,*)
         Write (*, '("Error(genfftmap): invalid ngrid : ", 3I8)') fftmap%ngrid
         Write (*,*)
         Stop
      End If
! total number of points in grid
      fftmap%ngrtot = fftmap%ngrid (1) * fftmap%ngrid (2) * fftmap%ngrid (3)
! determine integer ranges for grid
      intgv_custom (:, 1) = fftmap%ngrid (:) / 2 - fftmap%ngrid (:) + 1
      intgv_custom (:, 2) = fftmap%ngrid (:) / 2

      allocate(fftmap%igfft(ngrtot))

      Do ig = 1, ngrtot
         i1 = ivg (1, ig)
         i2 = ivg (2, ig)
         i3 = ivg (3, ig)
         
         If (i1 .Ge. 0) Then
            j1 = i1
         Else
            j1 = fftmap%ngrid (1) + i1
         End If
         
         If (i2 .Ge. 0) Then
            j2 = i2
         Else
            j2 = fftmap%ngrid (2) + i2
         End If
         If (i3 .Ge. 0) Then
            j3 = i3
         Else
            j3 = fftmap%ngrid (3) + i3
         End If
         if ((j1.ge.0).and.(j1.le.fftmap%ngrid (1)).and.(j2.ge.0).and.(j2.le.fftmap%ngrid (2)).and.(j3.ge.0).and.(j3.le.fftmap%ngrid (3))) then
           fftmap%igfft (ig) = j3 * fftmap%ngrid (2) * fftmap%ngrid (1) + j2 * fftmap%ngrid (1) + j1 + 1
         else
           fftmap%igfft (ig)=fftmap%ngrtot +1
         endif
      End Do

      Return
End Subroutine
!EOC
