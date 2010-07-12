!
! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
!
!
! !INTERFACE:
Subroutine gensymcmut (eps, maxsymcrys, nsymcrys, symlat, lsplsymc, &
& vtlsymc, scmut, tabel, tspainvsym)
! !DESCRIPTION:
!   Sets up the group multiplication table. The table is checked for consistency
!   in a way that it is required that every elements occurrs once and only once
!   in each row and column of the table. The first row and colmuns must consist
!   of the indentity since the first symmetry element is the identity by
!   convention.
!
! !REVISION HISTORY:
!   Created July 2008 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! arguments
      Real (8), Intent (In) :: eps
      Integer, Intent (In) :: maxsymcrys, nsymcrys
      Integer, Intent (In) :: symlat (3, 3, 48)
      Integer, Intent (In) :: lsplsymc (nsymcrys)
      Real (8), Intent (In) :: vtlsymc (3, maxsymcrys)
      Integer, Intent (Out) :: scmut (nsymcrys, nsymcrys)
      Logical, Intent (Out) :: tabel
      Logical, Intent (Out) :: tspainvsym
  ! local variables
      Integer :: isym, jsym, asym, lspli, lsplj, lspla, iv (3), nsis
      Integer :: doner (maxsymcrys), donec (maxsymcrys)
      Real (8) :: c (3, 3), ct (3, 3), si (3, 3), sj (3, 3), sa (3, 3), &
     & vtt (3), vtl (3), vtla (3)
      scmut (:, :) = 0
      nsis=0
      Do isym = 1, nsymcrys
         lspli = lsplsymc (isym)
         si (:, :) = dble (symlat(:, :, lspli))
         ! check for spatial inversion symmetry
         sj(:,:)=si(:,:) - &
         reshape((/-1.d0,0.d0,0.d0, 0.d0,-1.d0,0.d0, 0.d0,0.d0,-1.d0/), (/3,3/))
         if ((sum(abs(sj)) .lt. eps)) nsis=nsis+1
         Do jsym = 1, nsymcrys
            lsplj = lsplsymc (jsym)
            sj (:, :) = dble (symlat(:, :, lsplj))
        ! calculate rotation
            c = matmul (si, sj)
        ! first translation
            vtl (:) = matmul (sj, vtlsymc(:, jsym))
        ! second translation
            vtl (:) = vtl (:) + vtlsymc (:, isym)
            Call r3frac (eps, vtl, iv)
            vtl (:) = matmul (si, vtl)
            Call r3frac (eps, vtl, iv)
            Do asym = 1, nsymcrys
               lspla = lsplsymc (asym)
               sa (:, :) = dble (symlat(:, :, lspla))
           ! third translation
               vtla (:) = matmul (sa, vtlsymc(:, asym))
               Call r3frac (eps, vtla, iv)
           ! differece in rotation
               ct (:, :) = c (:, :) - sa (:, :)
           ! difference in translation
               vtt (:) = vtl (:) - vtla (:)
               If ((sum(Abs(ct)) .Lt. eps) .And. (sum(Abs(vtt)) .Lt. &
              & eps)) Then
              ! add element to multiplication table
                  scmut (isym, jsym) = asym
                  Cycle
               End If
            End Do
         End Do
      End Do
  ! check multiplication table for consistency
      Do isym = 1, nsymcrys
         donec (:) = 0
         doner (:) = 0
         Do jsym = 1, nsymcrys
            doner (scmut(isym, jsym)) = doner (scmut(isym, jsym)) + 1
            donec (scmut(jsym, isym)) = donec (scmut(jsym, isym)) + 1
         End Do
         Do jsym = 1, nsymcrys
            If (doner(jsym) .Ne. 1) Then
               Write (*,*)
               Write (*, '("Error(gensymcmut): error in multiplication &
              &table in row")')
               Write (*, '(" row number    : ", i6)') isym
               Write (*, '(" column number : ", i6)') jsym
               Write (*, '(" multiple occurrence : ", i6)') doner &
              & (jsym)
               Write (*,*)
               Stop
            End If
            If (donec(jsym) .Ne. 1) Then
               Write (*,*)
               Write (*, '("Error(gensymcmut): error in multiplication &
              &table in column")')
               Write (*, '(" row number    : ", i6)') jsym
               Write (*, '(" column number : ", i6)') isym
               Write (*, '(" multiple occurrence : ", i6)') donec &
              & (jsym)
               Write (*,*)
               Stop
            End If
         End Do
     ! check first row and column
         If ((scmut(1, isym) .Ne. isym) .Or. (scmut(isym, 1) .Ne. &
        & isym)) Then
            Write (*,*)
            Write (*, '("Error(gensymcmut): error in multiplication tab&
           &le")')
            Write (*, '(" first row or column have wrong property")')
            Write (*, '(" position : ", i6)') isym
            Write (*,*)
            Stop
         End If
      End Do
  ! check if group is Abelian
      tabel = .False.
      If (all((scmut-transpose(scmut)) .Eq. 0)) tabel = .True.
      tspainvsym=.false.
      ! check for spatial inversion symmetry
      if (nsis .gt. 0) tspainvsym=.true.
End Subroutine gensymcmut
!EOC
