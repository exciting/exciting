!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
Module summations
      Implicit None
Contains
!
!
      Subroutine doublesummation_simple_cz (z, za, zm, zb, alpha, beta, &
     & tbal)
         Implicit None
         Complex (8), Intent (Inout) :: z (:, :)
         Complex (8), Intent (In) :: za (:, :), zm (:, :), zb (:, :)
         Complex (8), Intent (In) :: alpha, beta
         Logical, Intent (In) :: tbal
    ! local variables
         Integer :: n1, n2, m1, m2, r1, r2, s1, s2
         Complex (8), Parameter :: zzero = cmplx (0.d0, 0.d0, 8), zone &
        & = cmplx (1.d0, 0.d0, 8)
         Complex (8), Allocatable :: zw (:, :)
         r1 = size (z, 1)
         r2 = size (z, 2)
         m1 = size (za, 1)
         n1 = size (za, 2)
         m2 = size (zb, 1)
         n2 = size (zb, 2)
         s1 = size (zm, 1)
         s2 = size (zm, 2)
    ! check output-array consistency
         If ((r1 .Ne. n1) .Or. (r2 .Ne. n2)) Then
            Write (*,*)
            Write (*, '("Error(doublesummation_simple_cz): inconsistent&
           & output array size: (r1,r2,n1,n2)",4i8)') r1, r2, n1, n2
            Write (*,*)
            Stop
         End If
    ! multiplication consistency
         If ((s1 .Ne. m1) .Or. (s2 .Ne. m2)) Then
            Write (*,*)
            Write (*, '("Error(doublesummation_simple_cz): inconsistent&
           & multiplication array size: (s1,s2,m1,m2)",4i8)') s1, s2, &
           & m1, m2
            Write (*,*)
            Stop
         End If
    ! choose the order of the two matrix multiplications in a way such that the
    ! number of operations is balanced for both multiplications
         If ((n1 .Ge. n2) .And. tbal) Then
            Allocate (zw(m1, n2))
       ! calculate double summation :
       ! Z(n1,n2) = alpha * sum{m1,m2} x
       !          x  conjg(A(m1,n1)) [M(m1,m2) B(m2,n2)] + beta*Z(n1,n2)
       ! as two matrix multiplications
       ! (i)  W(m1,n2) = M(m1,m2) . B(m2,n2)
       ! (ii) Z(n1,n2)= alpha * conjg(A'(n1,m1)) . W(m1,n2) + beta*Z(n1,n2)
       ! ad (i)
            Call zgemm ('n', 'n', m1, n2, m2, zone, zm, m1, zb, m2, &
           & zzero, zw, m1)
       ! ad (ii)
            Call zgemm ('c', 'n', n1, n2, m1, alpha, za, m1, zw, m1, &
           & beta, z, n1)
         Else
            Allocate (zw(n1, m2))
       ! calculate double summation :
       ! Z(n1,n2) = alpha * sum{m1,m2} x
       !          x  [conjg(A(m1,n1)) M(m1,m2)] B(m2,n2) + beta*Z(n1,n2)
       ! as two matrix multiplications
       ! (i)  W(n1,m2) = conjg(A'(n1,m1)) . M(m1,m2)
       ! (ii) Z(n1,n2)= alpha * W(n1,m2) .  B(m2,n2) + beta*Z(n1,n2)
       ! ad (i)
            Call zgemm ('c', 'n', n1, m2, m1, zone, za, m1, zm, m1, &
           & zzero, zw, n1)
       ! ad (ii)
            Call zgemm ('n', 'n', n1, n2, m2, alpha, zw, n1, zb, m2, &
           & beta, z, n1)
         End If
         Deallocate (zw)
      End Subroutine doublesummation_simple_cz
!
End Module summations
!
Module blaswrappers
      Implicit None
Contains
!
!
      Subroutine zgemm_wrap (z, ta, za, tb, zb, alpha, beta)
         Implicit None
    ! arguments
         Complex (8), Intent (Inout) :: z (:, :)
         Complex (8), Intent (In) :: za (:, :), zb (:, :)
         Character (1), Intent (In) :: ta, tb
         Complex (8), Intent (In) :: alpha, beta
    ! local variables
         Integer :: nz1, nz2, na1, na2, nb1, nb2
         nz1 = size (z, 1)
         nz2 = size (z, 2)
         If ((ta .Eq. 'n') .And. (tb .Eq. 'n')) Then
            na1 = size (za, 1)
            na2 = size (za, 2)
            nb1 = size (zb, 1)
            nb2 = size (zb, 2)
       ! z(nz1,nz2)~alpha*za(na1,na2)*zb(nb1,nb2)
            If ((na1 .Gt. nz1) .Or. (na2 .Ne. nb1) .Or. (nb2 .Gt. nz2)) &
           & Then
               Write (*,*)
               Write (*, '("Error(zgemm_wrap): inconsistent array dimen&
              &sions")')
               Write (*, '(" Z = alpha A B + beta Z")')
               Write (*, '(" Z : ",2i8)') nz1, nz2
               Write (*, '(" A : ",2i8)') na1, nb2
               Write (*, '(" B : ",2i8)') nb1, nb2
               Write (*,*)
               Stop
            End If
       ! call to BLAS routine
            Call zgemm (ta, tb, na1, nb2, na2, alpha, za, na1, zb, nb1, &
           & beta, z, nz1)
         Else If ((ta .Eq. 'n') .And. (tb .Ne. 'n')) Then
            na1 = size (za, 1)
            na2 = size (za, 2)
            nb1 = size (zb, 1)
            nb2 = size (zb, 2)
       ! z(nz1,nz2)~alpha*za(na1,na2)* j(zb(nb2,nb1))
            If ((na1 .Gt. nz1) .Or. (na2 .Ne. nb2) .Or. (nb1 .Gt. nz2)) &
           & Then
               Write (*,*)
               Write (*, '("Error(zgemm_wrap): inconsistent array dimen&
              &sions")')
               Write (*, '(" Z = alpha A j(B) + beta Z")')
               Write (*, '(" Z    : ",2i8)') nz1, nz2
               Write (*, '(" A    : ",2i8)') na1, nb2
               Write (*, '(" j(B) : ",2i8)') nb2, nb1
               Write (*,*)
               Stop
            End If
       ! call to BLAS routine
            Call zgemm (ta, tb, na1, nb1, na2, alpha, za, na1, zb, nb1, &
           & beta, z, nz1)
         Else If ((ta .Ne. 'n') .And. (tb .Eq. 'n')) Then
            na1 = size (za, 1)
            na2 = size (za, 2)
            nb1 = size (zb, 1)
            nb2 = size (zb, 2)
       ! z(nz1,nz2)~alpha*j(za(na2,na1))*zb(nb1,nb2)
            If ((na2 .Gt. nz1) .Or. (na1 .Ne. nb1) .Or. (nb2 .Gt. nz2)) &
           & Then
               Write (*,*)
               Write (*, '("Error(zgemm_wrap): inconsistent array dimen&
              &sions")')
               Write (*, '(" Z = alpha A B + beta Z")')
               Write (*, '(" Z    : ",2i8)') nz1, nz2
               Write (*, '(" j(A) : ",2i8)') na2, nb1
               Write (*, '(" B    : ",2i8)') nb1, nb2
               Write (*,*)
               Stop
            End If
       ! call to BLAS routine
            Call zgemm (ta, tb, na2, nb2, na1, alpha, za, na1, zb, nb1, &
           & beta, z, nz1)
         Else If ((ta .Ne. 'n') .And. (tb .Ne. 'n')) Then
            na1 = size (za, 1)
            na2 = size (za, 2)
            nb1 = size (zb, 1)
            nb2 = size (zb, 2)
       ! z(nz1,nz2)~alpha*j(za(na2,na1)) * j(zb(nb2,nb1))
            If ((na2 .Gt. nz1) .Or. (na1 .Ne. nb2) .Or. (nb1 .Gt. nz2)) &
           & Then
               Write (*,*)
               Write (*, '("Error(zgemm_wrap): inconsistent array dimen&
              &sions")')
               Write (*, '(" Z = alpha A j(B) + beta Z")')
               Write (*, '(" Z    : ",2i8)') nz1, nz2
               Write (*, '(" j(A) : ",2i8)') na2, nb1
               Write (*, '(" j(B) : ",2i8)') nb2, nb1
               Write (*,*)
               Stop
            End If
       ! call to BLAS routine
            Call zgemm (ta, tb, na2, nb1, na1, alpha, za, na1, zb, nb1, &
           & beta, z, nz1)
         End If
      End Subroutine zgemm_wrap
End Module blaswrappers
