!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: genwiqggp
! !INTERFACE:
!
!
Subroutine genwiqggp (flag, iq, igq1, igq2, clwt)
! !USES:
      Use modmain
      Use modxs
      Use m_genfilname
      Use m_getunit
! !DESCRIPTION:
!   Effective integrals of Coulomb interaction. See routine {\tt genwiq2}.
!
! !REVISION HISTORY:
!   Created February 2008 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! arguments
      Integer, Intent (In) :: flag, iq, igq1, igq2
      Real (8), Intent (Out) :: clwt
  ! local variables
      Integer, Parameter :: ns0 = 10, nss = 20
      Integer :: np, ns, i1, i2, i3, i, ip, nrbox
      Real (8) :: d (3), dv, sum2, t1, t11, t22, t33, t2
      Real (8) :: blim (2), blen, vllim (3), ran (3), ranl (3), &
     & omegabox
      Real (8) :: qsz
      Real (8), Allocatable :: xa (:), ya (:), c (:)
      Real (8) :: vsc (3), v01 (3), v02 (3), v11 (3), v21 (3), v31 (3), &
     & v32 (3)
  ! external functions
      Real (8) :: polynom
      External polynom
  ! determine G+q-vectors
      v01 (:) = vgqc (:, igq1, iq)
      v02 (:) = vgqc (:, igq2, iq)
  ! integrate out singularites and/or improve accuracy in summations
      Select Case (flag)
      Case (0)
     ! modified Coulomb potential in reciprocal space
         clwt = sptclg (igq1, iq) * sptclg (igq2, iq)
      Case (1)
     ! radius of sphere with same volume than subcell
         qsz = (6*pi**2/(omega*nqpt)) ** (1.d0/3.d0)
         If ((igq1 .Ne. 1) .Or. (igq2 .Ne. 1)) Then
        ! integrate out 1/q singularity by spherical Volume
            clwt = (qsz**2*omega*nqpt/pi) / gqc (igq2, iq)
         Else If ((igq1 .Eq. 1) .And. (igq2 .Eq. 1)) Then
        ! integrate out 1/q^2 singularity by spherical Volume
            clwt = 2 * qsz * omega * nqpt / pi
         Else
            Write (*,*)
            Write (*, '("Error(genwiqggp): analytic method chosen for r&
           &egular case")')
            Write (*,*)
            Call terminate
         End If
      Case (2)
         np = 2
     ! higher order extrapolation for 1/q^2 term
         If ((igq1 .Eq. 1) .And. (igq2 .Eq. 1)) np = 3
         Allocate (xa(np), ya(np), c(np))
     ! loop over different subdivisions
         ns = ns0
         Do ip = 1, np
        ! subdivision vectors in lattice coordinates
            Do i = 1, 3
               d (i) = 1.d0 / (dble(ngridq(i)*2*ns))
            End Do
        ! smallest volume element (we drop the (2pi)^3/omega factor!)
            dv = d (1) * d (2) * d (3)
        ! compute the integral of 1/(|p+q+G||p+q+Gp|) over the small
        ! fraction of the Brillouin zone centered at the Gamma-point
            sum2 = 0.d0
        ! the p-point grid is started here
            Do i1 = - ns, ns - 1
               t11 = dble (i1) * d (1)
               Do i2 = - ns, ns - 1
                  t22 = dble (i2) * d (2)
                  Do i3 = - ns, ns - 1
                     t33 = dble (i3) * d (3)
                 ! p-vector
                     vsc (:) = t11 * bvec (:, 1) + t22 * bvec (:, 2) + &
                    & t33 * bvec (:, 3)
                 ! p+q+G vector
                     v31 (:) = v01 (:) + vsc (:)
                 ! p+q+Gp vector
                     v32 (:) = v02 (:) + vsc (:)
                     t2 = Sqrt (sum(v31**2)*sum(v32**2))
                 ! check if integrand would approach singularity
                     If (t2 .Gt. 1.d-14) Then
                        sum2 = sum2 + 1.d0 / t2
                     End If
                  End Do
               End Do
            End Do
            sum2 = sum2 * dv
            xa (ip) = dv ** (1.d0/3.d0)
            ya (ip) = sum2
        ! increment number of subdivisions
            ns = ns + nss
         End Do
     ! extrapolate the volume element to zero with a polynomial
         clwt = polynom (0, np, xa, ya, c, 0.d0) * fourpi * nqpt
         Deallocate (xa, ya, c)
      Case (3)
     ! RIM method: documentation of Self code by Andrea Marini
     ! find maximum extension of small Brillouin zone
         blim (:) = 0.d0
         Do i1 = - 1, 1, 2
            v11 (:) = dble (i1) * bvec (:, 1) / dble (ngridq(1)*2)
            Do i2 = - 1, 1, 2
               v21 (:) = v11 (:) + dble (i2) * bvec (:, 2) / dble &
              & (ngridq(2)*2)
               Do i3 = - 1, 1, 2
                  v31 (:) = v21 (:) + dble (i3) * bvec (:, 3) / dble &
                 & (ngridq(3)*2)
                  Do i = 1, 3
	         ! lower limit for box
                     If (v31(i) .Lt. blim(1)) blim (1) = v31 (i)
		 ! upper limit for box
                     If (v31(i) .Gt. blim(2)) blim (2) = v31 (i)
                  End Do
               End Do
            End Do
         End Do
     ! box length
         blen = blim (2) - blim (1)
     ! limits of sBZ in lattice
         vllim (:) = 1.d0 / dble (2*ngridq(:))
     ! needs high value (above 1e6) to converge - inferior to method Nr. 2
         nrbox = 1000000
         t1 = 0.d0
         omegabox = blen ** 3
         Do i = 1, 10000
            Call random_number (t2)
         End Do
         Do i = 1, nrbox
            Call random_number (ran)
	! map random number to box
            ran (:) = ran (:) * blen
            ran (:) = ran (:) + blim (1)
	! check if random vector is in sBZ
            ranl = matmul (binv, ran)
            If (all(ranl .Gt.-vllim) .And. (all(ranl .Lt. vllim)) .And. &
           & (sum(Abs(ran)) .Gt. 1.d-14)) Then
               t2 = sum ((v01+ran)**2) * sum ((v02+ran)**2)
               t1 = t1 + 1.d0 / Sqrt (t2)
            End If
         End Do
         clwt = t1 * (omegabox/nrbox) * fourpi * nqpt * omega / (twopi) &
        & ** 3
      End Select
End Subroutine genwiqggp
!EOC
