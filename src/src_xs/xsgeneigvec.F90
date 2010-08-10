!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine xsgeneigvec
      Use modmain
      Use modinput
      Use modmpi
      Use modxs
      Use m_writegqpts
      Use m_filedel
      Use m_genfilname
      Implicit None
  ! local variables
      Character (*), Parameter :: thisnam = 'xsgeneigvec'
      Real (8) :: vqlt (3)
      Integer :: iq, qi, qf
      Logical, External :: tqgamma
  ! initialize universal variables
      Call init0
      Call init1
      Call init2
  ! SCF allready parallelized for k-point set
  ! add Gamma Q-point
      qi = 0
      qf = nqpt
  ! if first Q-point is Gamma-point we copy files
      If (tqgamma(1)) qi = 1
  ! for screening there is just the Gamma q-point
      If (tscreen) Then
         qi = 0
         qf = 0
      End If
  ! write q-points
      If (rank .Eq. 0) Call writeqpts
  ! calculate eigenvectors for each q-point (k+q point set)
      Do iq = qi, qf
         If (.Not. tscreen) Call genfilname (iqmt=Max(0, iq), setfilext=.True.)
         vqlt (:) = 0.d0
         If (iq .Ne. 0) then
            vqlt (:) = vql (:, iq)
            If (rank .Eq. 0) Call writegqpts (iq, filext)
         End If
     ! write eigenvectors, -values, occupancies and contracted MT coefficients
         Call writeevec (vqlt, qvkloff(1, iq), filext)
         If (.Not. tscreen) Then
            Write (unitout, '("Info(", a, "): eigenvectors generated fo&
           &r Q-point (iq, vql below)")') thisnam
            Write (unitout, '(i6, 3g18.10)') iq, vqlt (:)
         Else
            Write (unitout, '(a)') 'Info(' // thisnam // '): eigenvecto&
           &rs generated for associated(input%xs%screening)/screened in&
           &teraction:'
         End If
         If (rank .Eq. 0) Then
        ! safely remove unnecessary files
            Call filedel ('EQATOMS'//trim(filext))
            Call filedel ('EVALCORE'//trim(filext))
            Call filedel ('FERMIDOS'//trim(filext))
            Call filedel ('GEOMETRY'//trim(filext))
            Call filedel ('LATTICE'//trim(filext))
            Call filedel ('IADIST'//trim(filext))
            Call filedel ('LINENGY'//trim(filext))
            Call filedel ('SYMCRYS'//trim(filext))
            Call filedel ('SYMLAT'//trim(filext))
            Call filedel ('SYMSITE'//trim(filext))
            Call filedel ('TOTENERGY'//trim(filext))
            Call filedel ('EVALFV'//trim(filext))
            Call filedel ('RMSDVEFF'//trim(filext))
            Call filedel ('DTOTENERGY'//trim(filext))
            if (input%groundstate%tforce) Call filedel ('DFORCEMAX'//trim(filext))
            Call filedel ('CHGDIST'//trim(filext))
            Call filedel ('SYMGENR'//trim(filext))
            Call filedel ('SYMINV'//trim(filext))
            Call filedel ('SYMMULT'//trim(filext))
            Call filedel ('SYMMULT_TABLE'//trim(filext))
            Call filedel ('SYMT2'//trim(filext))
         End If
     ! end loop over q-points
      End Do
      If ((rank .Eq. 0) .And. tqgamma(1) .And. ( .Not. tscreen)) Then
         Write (unitout, '("Info(", a, "): first Q-point is Gamma-point&
        & - copying relevant files")') thisnam
     ! write files again one by one
         Call copyfilesq0
      End If
      Call barrier
      Write (unitout, '("Info(", a, "): generation of eigenvectors fini&
     &shed")') thisnam
      Call genfilname (setfilext=.True.)
End Subroutine xsgeneigvec
