!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine findocclims (iq, iocc0, iocc, iunocc0, iunocc, io0, io, iu0, &
& iu)
      Use modmain
      Use modinput
      Use modxs
      Use m_genfilname
      Implicit None
  ! arguments
      Integer, Intent (In) :: iq
      Integer, Intent (Out) :: iocc0, iocc, iunocc0, iunocc
      Integer, Intent (Out) :: io0 (nkpt), io (nkpt), iu0 (nkpt), iu &
     & (nkpt)
  ! local variables
      Integer :: ik, ikq, i0, i
      Logical :: t
      t = allocated (evalsv0)
      If ( .Not. t) allocate (evalsv0(nstsv, nkpt))
      Do ik = 1, nkpt
     ! k+q-point set
         ikq = ik
         If (iq .Ne. 0) ikq = ikmapikq (ik, iq)
         Call getoccsv (vkl(1, ikq), occsv(1, ikq))
         Call getevalsv (vkl(1, ikq), evalsv(1, ikq))
         Do i = 1, nstsv
            If (occsv(i, ikq) .Lt. input%groundstate%epsocc) Exit
         End Do
         io (ik) = i - 1
         Do i = nstsv, 1, - 1
            If (occsv(i, ikq) .Gt. (occmax-input%groundstate%epsocc)) &
           & Exit
         End Do
         iu (ik) = i + 1
         If (iq .Ne. 0) Then
        ! k-point set (q=0)
            Call getoccsv0 (vkl0(1, ik), occsv0(1, ik))
            Call getevalsv0 (vkl0(1, ik), evalsv0(1, ik))
            Do i0 = 1, nstsv
               If (occsv0(i0, ik) .Lt. input%groundstate%epsocc) Exit
            End Do
            io0 (ik) = i0 - 1
            Do i0 = nstsv, 1, - 1
               If (occsv0(i0, ik) .Gt. &
              & (occmax-input%groundstate%epsocc)) Exit
            End Do
            iu0 (ik) = i0 + 1
         Else
            io0 (ik) = io (ik)
            iu0 (ik) = iu (ik)
         End If
      End Do
      If (iq .Ne. 0) Then
     ! lowest and highest valence energy
         evlmin = Min (minval(evalsv(1, :)), minval(evalsv0(1, :)))
         evlmax = Max (maxval(evalsv(nstsv, :)), maxval(evalsv0(nstsv, &
        & :)))
     ! lower and higher cutoff valence energy
         evlmincut = Max (maxval(evalsv(1, :)), maxval(evalsv0(1, :)))
         evlmaxcut = Min (minval(evalsv(nstsv, :)), &
        & minval(evalsv0(nstsv, :)))
      Else
     ! lowest and highest valence energy
         evlmin = minval (evalsv(1, :))
         evlmax = maxval (evalsv(nstsv, :))
     ! lower and higher cutoff valence energy
         evlmincut = maxval (evalsv(1, :))
         evlmaxcut = minval (evalsv(nstsv, :))
      End If
  ! overall highest (partially) occupied state
      iocc0 = maxval (io0)
      iocc = maxval (io)
  ! overall lowest (partially) unoccupied state
      iunocc0 = minval (iu0)
      iunocc = minval (iu)
  ! the maximum/minimum value is used since a shifted (k+q)-mesh which is not
  ! commensurate can cause partially occupied states that are absent for the
  ! k-mesh
      iocc0 = Max (iocc0, iocc)
      iocc = iocc0
      iunocc0 = Min (iunocc0, iunocc)
      iunocc = iunocc0
  ! determine if system has a gap in energy
      If (iq .Ne. 0) Then
     ! highest (partially) occupied state energy
         evlhpo = Max (maxval(evalsv(iocc0, :)), maxval(evalsv0(iocc0, &
        & :)))
     ! lowest (partially) unoccupied state energy
         evllpu = Min (minval(evalsv(iunocc0, :)), &
        & minval(evalsv0(iunocc0, :)))
      Else
     ! highest (partially) occupied state energy
         evlhpo = maxval (evalsv(iocc0, :))
     ! lowest (partially) unoccupied state energy
         evllpu = minval (evalsv(iunocc0, :))
      End If
  ! determine if system has a gap in energy
      ksgap = evlhpo .Lt. efermi
  ! assign nstocc0 and nstunocc0
      nstocc0 = iocc0
      nstunocc0 = nstsv - nstocc0
      If ((iocc0 .Ge. iunocc) .Or. (iocc .Ge. iunocc0)) Then
         Write (unitout, '(a)') 'Info(findocclims): partially occupied &
        &states present'
      End If
      If (ksgap) Then
         Write (unitout, '(a)') 'Info(findocclims): system has Kohn-Sha&
        &m gap'
      Else
         Write (unitout, '(a)') 'Info(findocclims): no Kohn-Sham gap fo&
        &und'
      End If
  ! debug output
      If (input%xs%dbglev .Gt. 0) Then
         Write (*, '(a)') 'Debug(findocclims):'
         Write (*, '(a)') ' iocc0, iocc, iunocc0, iunocc below:'
         Write (*, '(4i8)') iocc0, iocc, iunocc0, iunocc
         Write (*, '(a)') ' ik, io0, iu, diff, io, iu0, diff below:'
         Do ik = 1, nkpt
            Write (*, '(7i8)') ik, io0 (ik), iu (ik), iu (ik) - io0 &
           & (ik), io (ik), iu0 (ik), iu0 (ik) - io (ik)
         End Do
         Write (*,*)
      End If
      If ( .Not. t) deallocate (evalsv0)
End Subroutine findocclims
