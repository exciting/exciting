
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine mixerifc (mtype, n, v, dv, mode)
!use modmain
      Use modinput
      Use mod_potential_and_density, Only:
      Use mod_spin, Only:
      Use mod_convergence, Only: iscl
      Use modmixermsec
      Use modmixadapt
      Implicit None
! arguments
      Integer, Intent (In) :: mtype
      Integer, Intent (In) :: n
      Real (8), Intent (Inout) :: v (n)
      Real (8), Intent (Out) :: dv
      Integer, Intent (Inout) :: mode
!mode: 	-1 call initialisation routines,
!		-2:call destructor
!		else ignore
      Integer, Parameter :: maxsd = 3
      Select Case (mtype)
      Case (1)
! adaptive linear mixing
! calculate memmory requirement if mode negative
         If (mode .Eq.-1) Then
            mode = 0
            If (allocated(work)) deallocate (work)
            Allocate (work(3*n))
!            Return
         End If
         If (mode .Eq.-2) Then
            Deallocate (work)
            Return
         End If
!--
         Call mixadapt (iscl, input%groundstate%beta0, &
        & input%groundstate%betainc, input%groundstate%betadec, n, v, &
        & work, work(n+1), work(2*n+1), dv)
!
      Case (2)
 ! multicecant broyden
         If (mode .Eq.-1) Then
            Call initmixermsec (n,input%groundstate%msecStoredSteps)
            mode = 0
!            Return
         End If
         If (mode .Eq.-2) Then
            Call freearraysmixermsec ()
            Return
         End If
         Call mixmsec (iscl, v, dv, n)
      Case (3)
 ! Pulay mixing
         If (associated(input%groundstate%spin)) Then
            call warning('Warning(mixerifc): Pulay mixing problematic with spin-polarised calculations')
         End If
         If (mode .Eq.-1) Then
            Allocate (work(2*n*maxsd))
            mode = 0
            Return
         End If
         If (mode .Eq.-2) Then
            Deallocate(work)
            Return
         End If
         Call mixpulay(iscl, n, maxsd, v, work, work(n*maxsd+1), dv)
      Case Default
         Write(*,*)
         Write(*, '("Error(mixerifc): mtype not defined : ", I8)') mtype
         Write(*,*)
         Stop
      End Select
      Return
End Subroutine
