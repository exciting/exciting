module m_wsweight
    implicit none
    contains
!
!BOP
! !ROUTINE: ws_weight
! !INTERFACE:
!
!
Subroutine ws_weight (vrl, vrsl, vpl, zwght, kgrid)
! !INPUT/OUTPUT PARAMETERS:
!   vrl   : R vector in lattice coordinates (in,real(3))
!   vrsl  : R+s vector (vector from atom 1 in cell at origin to atom 2 in
!           unit cell at R) in lattice coordinates (in,real(3))
!   vpl   : q vector in reciprocal lattice coordinates (in,real(3))
!   zwght : combined complex weight to be used in the Fourier 
!           transformation (out,complex)
! !DESCRIPTION:
!   Folds the ${\bf R}+{\bf s}$ vector back into the Wigner-Seitz cell of the supercell, 
!   ${\bf s}$ being the atom basis. \\
!   Then equivalent vectors (total number $n$) are searched for and the combined weights and 
!   phase factors for the Fourier transformation computed and returned:
!   $$ \frac{1}{n}\sum_{j=1}^{n} e^{-i {\bf q}\cdot {\bf R}_j} $$
!
! !REVISION HISTORY:
!   Created September 2012 (STK)
!EOP
!BOC
      Use modmain
      Implicit None
! arguments
      Real(8), Intent(In) :: vrl(3)
      Real(8), Intent(In) :: vrsl(3)
      Real(8), Intent(In) :: vpl(3)
      Complex(8), Intent(Out) :: zwght
      logical, intent( in), optional :: kgrid
! local variables
      Integer :: i1, i2, i3, j1, j2, j3, n, ngrid_(3)
      Real(8) :: v0(3), v1(3), v2(3), v3(3), t1, t2, t3
      Real(8) :: vrssl(3), vsl(3), vr2l(3), vrs2l(3)
      Integer :: ivrl(3)
      ngrid_ = ngridq
      if( present( kgrid)) then
        if( kgrid) ngrid_ = input%groundstate%ngridk
      end if
      vsl(:) = vrsl(:) - vrl(:)
! vrssl is vector to atom s in cell R; referring to supercell
      vrssl(:) = vrsl(:) / dble(ngrid_(:))
! map vector to [0,1) interval
      Call r3frac (input%structure%epslat, vrssl, ivrl)
! convert to cartesian and shift by lattice vectors to obtain shortest possible vector
      Call r3mv (input%structure%crystal%basevect, vrssl, v0)
      t1 = v0(1)**2 + v0(2)**2 + v0(3)**2
      Do i1 = -1,1
         v1(:) = v0(:) + dble(i1)*input%structure%crystal%basevect(:,1)
         Do i2 = -1,1
            v2(:) = v1(:) + dble(i2)*input%structure%crystal%basevect(:,2)
            Do i3 = -1,1
               v3(:) = v2(:) + dble(i3)*input%structure%crystal%basevect(:,3)
               t2 = v3(1)**2 + v3(2)**2 + v3(3)**2
               If (t2 .Lt. (t1 + input%structure%epslat)) Then
                  j1 = i1
                  j2 = i2
                  j3 = i3
                  t1 = t2
                  v0(:) = v3(:)
               End If
            End Do
         End Do
      End Do
! search equivalent vectors
      n = 0
      zwght = (0.d0, 0.d0)
      Do i1 = -1,1
         v1(:) = v0(:) + dble(i1)*input%structure%crystal%basevect(:,1)
         Do i2 = -1,1
            v2(:) = v1(:) + dble(i2)*input%structure%crystal%basevect(:,2)
            Do i3 = -1,1
               v3(:) = v2(:) + dble(i3)*input%structure%crystal%basevect(:,3)
               t2 = v3(1)**2 + v3(2)**2 + v3(3)**2
               If (abs(t2 - t1) .Lt. input%structure%epslat) Then
                  n = n + 1
! convert to lattice coords of original unit cell and consider R part only
                  Call r3mv (ainv, v3, vrssl)
                  vrs2l(:) = vrssl(:)*dble(ngrid_(:))
                  vr2l(:) = vrs2l(:) - vsl(:)
                  t3 = -twopi*(vpl(1)*vr2l(1) + vpl(2)*vr2l(2) &
                           & + vpl(3)*vr2l(3))
! sum up phase factors
                  zwght = zwght + cmplx(cos(t3), sin(t3), 8)
               End If
            End Do
         End Do
      End Do
      call r3mv( ainv, v0, v1)
! divide by number of equivalent vectors
      zwght = zwght / dble(n)
      Return
End Subroutine ws_weight
!EOC
end module m_wsweight
