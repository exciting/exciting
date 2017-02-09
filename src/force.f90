!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: force
! !INTERFACE:
!
!
Subroutine force
! !USES:
      Use modinput
      Use modmain
      Use modmpi
! !DESCRIPTION:
!   Computes the various contributions to the atomic forces. In principle, the
!   force acting on a nucleus is simply the gradient at that site of the
!   classical electrostatic potential from the other nuclei and the electronic
!   density. This is a result of the Hellmann-Feynman theorem. However because
!   the basis set is dependent on the nuclear coordinates and is not complete,
!   the Hellman-Feynman force is inacurate and corrections to it are required.
!   The first is the core correction which arises because the core wavefunctions
!   were determined by neglecting the non-spherical parts of the effective
!   potential $v_s$. Explicitly this is given by
!   $$ {\bf F}_{\rm core}^{\alpha}=\int_{\rm MT_{\alpha}} v_{\rm s}({\bf r})
!    \nabla\rho_{\rm core}^{\alpha}({\bf r})\,d{\bf r} $$
!   for atom $\alpha$. The second, which is the incomplete basis set (IBS)
!   correction, is due to the position dependence of the APW functions, and is
!   derived by considering the change in total energy if the eigenvector
!   coefficients were fixed and the APW functions themselves were changed. This
!   would result in changes to the first-variational Hamiltonian and overlap
!   matrices given by
!   \begin{align*}
!    \delta H_{\bf G,G'}^{\alpha}&=i({\bf G-G'})
!    \left(H^{\alpha}_{\bf G+k,G'+k}-\frac{1}{2}({\bf G+k})\cdot({\bf G'+k})
!    \tilde{\Theta}_{\alpha}({\bf G-G'})e^{-i({\bf G-G'})\cdot{\bf r}_{\alpha}}
!    \right)\\
!    \delta O_{\bf G,G'}^{\alpha}&=i({\bf G-G'})\left(O^{\alpha}_{\bf G+k,G'+k}
!    -\tilde{\Theta}_{\alpha}({\bf G-G'})e^{-i({\bf G-G'})\cdot{\bf r}_{\alpha}}
!    \right)
!   \end{align*}
!   where both ${\bf G}$ and ${\bf G'}$ run over the APW indices;
!   $\tilde{\Theta}_{\alpha}$ is the form factor of the smooth step function for
!   muffin-tin $\alpha$; and $H^{\alpha}$ and $O^{\alpha}$ are the muffin-tin
!   Hamiltonian and overlap matrices, respectively. The APW-local-orbital part
!   is given by
!   \begin{align*}
!    \delta H_{\bf G,G'}^{\alpha}&=i({\bf G+k})H^{\alpha}_{\bf G+k,G'+k}\\
!    \delta O_{\bf G,G'}^{\alpha}&=i({\bf G+k})O^{\alpha}_{\bf G+k,G'+k}
!   \end{align*}
!   where ${\bf G}$ runs over the APW indices and ${\bf G'}$ runs over the
!   local-orbital indices. There is no contribution from the
!   local-orbital-local-orbital part of the matrices. We can now write the IBS
!   correction in terms of the basis of first-variational states as
!   \begin{align*}
!    {\bf F}_{ij}^{\alpha{\bf k}}=\sum_{\bf G,G'}
!    b^{i{\bf k}*}_{\bf G}b^{j{\bf k}}_{\bf G'}\left(
!    \delta H_{\bf G,G'}^{\alpha}-\epsilon_j\delta O_{\bf G,G'}^{\alpha}\right),
!   \end{align*}
!   where $b^{i{\bf k}}$ is the first-variational eigenvector.
!   Finally, the ${\bf F}_{ij}^{\alpha{\bf k}}$ matrix elements can be
!   multiplied by the second-variational coefficients, and contracted over all
!   indices to obtain the IBS force:
!   \begin{align*}
!    {\bf F}_{\rm IBS}^{\alpha}=\sum_{\bf k}w_{\bf k}\sum_{l\sigma}n_{l{\bf k}}
!    \sum_{ij}c_{\sigma i}^{l{\bf k}*}c_{\sigma j}^{l{\bf k}}
!    {\bf F}_{ij}^{\alpha{\bf k}}
!    +\int_{\rm MT_{\alpha}}v_{\rm s}({\bf r})\nabla\left[\rho({\bf r})
!    -\rho^{\alpha}_{\rm core}({\bf r})\right]\,d{\bf r},
!   \end{align*}
!   where $c^{l{\bf k}}$ are the second-variational coefficients, $w_{\bf k}$
!   are the $k$-point weights, $n_{l{\bf k}}$ are the occupancies, and
!   $v_{\rm s}$ is the Kohn-Sham effective potential. See routines {\tt hmlaa},
!   {\tt olpaa}, {\tt hmlalo}, {\tt olpalo}, {\tt energy}, {\tt seceqn} and
!   {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!   Fixed problem with second-variational forces, May 2008 (JKD)
!   k-point parallelisation of IBS forces, October 2013 (Andris)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: ik, is, ia, ias, nr, i
      Real (8) :: sum, t1
      Real (8) :: ts0, ts1
! allocatable arrays
      Real (8), Allocatable :: rfmt (:, :)
      Real (8), Allocatable :: grfmt (:, :, :)
      Real (8), Allocatable :: ffacg (:, :)
      Real (8), Allocatable :: forcesum(:,:)
! external functions
      Real (8) :: rfmtinp
      External rfmtinp
      Call timesec (ts0)

      If (allocated(forcetot)) deallocate (forcetot)
      Allocate (forcetot(3, natmtot))

      Allocate (rfmt(lmmaxvr, nrmtmax))
      Allocate (grfmt(lmmaxvr, nrmtmax, 3))
!--------------------------------!
!     Hellmann-Feynman force     !
!--------------------------------!
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
! compute the gradient of the Coulomb potential
            Call gradrfmt (1, nrmt(is), spr(:, is), lmmaxvr, nrmtmax, &
           & vclmt(:, :, ias), grfmt)
            forcehf (:, ias) = - spzn (is) * grfmt (1, 1, :) * y00
         End Do
      End Do
! symmetrise Hellmann-Feynman force
      Call symvect (.False., forcehf)
!--------------------------------------!
!     core correction to the force     !
!--------------------------------------!
      rfmt (:, :) = 0.d0
      Do is = 1, nspecies
         nr = nrmt (is)
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
! compute the gradient of the core density
            rfmt (1, 1:nr) = rhocr (1:nr, ias) / y00
            Call gradrfmt (1, nr, spr(:, is), lmmaxvr, nrmtmax, rfmt, &
           & grfmt)
            Do i = 1, 3
               forcecr (i, ias) = rfmtinp (1, 1, nr, spr(:, is), &
              & lmmaxvr, veffmt(:, :, ias), grfmt(:, :, i))
            End Do
         End Do
      End Do
! symmetrise core correction force
      Call symvect (.False., forcecr)
!-------------------------------------!
!     IBS correction to the force     !
!-------------------------------------!
! set the IBS forces to zero
      forceibs (:, :) = 0.d0



      If (input%groundstate%tfibs) Then
         Allocate (ffacg(ngvec, nspecies))
! generate the step function form factors
         Do is = 1, nspecies
            Call genffacg (is, ffacg(:, is))
         End Do

         call olprad

! compute k-point dependent contribution to the IBS force

#ifdef MPI
         Do ik = firstk (rank), lastk (rank)
#else
         Do ik = 1, nkpt
#endif
            Call forcek (ik, ffacg)
         End Do

#ifdef MPI
        allocate(forcesum(3,natmtot))
        forcesum=0d0
        call MPI_ALLREDUCE(forceibs, forcesum, natmtot*3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        forceibs=forcesum
        deallocate(forcesum)
#endif

! integral of effective potential with gradient of valence density
         Do is = 1, nspecies
            nr = nrmt (is)
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               rfmt (:, 1:nr) = rhomt (:, 1:nr, ias)
               rfmt (1, 1:nr) = rfmt (1, 1:nr) - rhocr (1:nr, ias) / &
              & y00
               Call gradrfmt (input%groundstate%lmaxvr, nr, spr(:, is), &
              & lmmaxvr, nrmtmax, rfmt, grfmt)
               Do i = 1, 3
                  t1 = rfmtinp (1, input%groundstate%lmaxvr, nr, spr(:, &
                 & is), lmmaxvr, veffmt(:, :, ias), grfmt(:, :, i))
                  forceibs (i, ias) = forceibs (i, ias) + t1
               End Do
            End Do
         End Do

! symmetrise IBS force
         Call symvect (.False., forceibs)
         Deallocate (ffacg)
      End If
! remove net forces of ibs, core, and HF forces (center of mass should not move)

      Do i = 1, 3
         sum = 0.d0
         Do ias = 1, natmtot
            sum = sum + forceibs (i, ias)
         End Do
         sum = sum / dble (natmtot)
         forceibs (i, :) = forceibs (i, :) - sum
      End Do

      Do i = 1, 3
         sum = 0.d0
         Do ias = 1, natmtot
            sum = sum + forcecr (i, ias)
         End Do
         sum = sum / dble (natmtot)
         forcecr (i, :) = forcecr (i, :) - sum
      End Do

      Do i = 1, 3
         sum = 0.d0
         Do ias = 1, natmtot
            sum = sum + forcehf (i, ias)
         End Do
         sum = sum / dble (natmtot)
         forcehf (i, :) = forcehf (i, :) - sum
      End Do

! calculate total force

      Do ias = 1, natmtot
         forcetot (:, ias) = forcehf (:, ias) + forcecr (:, ias) + forceibs (:, ias)
      End Do

!dispersion correction

      If ( tlast .And. input%groundstate%vdWcorrection .Ne. "none" ) Then
         If (allocated(force_disp)) deallocate (force_disp)
         Allocate(force_disp(3,natmtot))
         If ( input%groundstate%vdWcorrection .Eq. "DFTD2" ) Then
            Call DFT_D2_force
         Else If ( input%groundstate%vdWcorrection .Eq. "TSvdW" ) Then
            Call TS_vdW_force
         End If
         Do ias = 1, natmtot
            forcetot (:, ias) = forcetot (:, ias) + force_disp (:, ias)
         End Do
      End If

! symmetrise total force

      Call symvect (.False., forcetot)

! remove net total force (center of mass should not move)

      Do i = 1, 3
         sum = 0.d0
         Do ias = 1, natmtot
            sum = sum + forcetot (i, ias)
         End Do
         sum = sum / dble (natmtot)
         forcetot (i, :) = forcetot (i, :) - sum
      End Do

! compute maximum force magnitude over all non constrained atoms and components

      forcemax = 0.d0
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            t1 = 0.d0
            Do i = 1, 3
               if (.not.input%structure%speciesarray(is)%species%atomarray(ia)%atom%lockxyz(i)) &
              &   t1 = t1 + forcetot(i, ias)**2
            End Do
            t1 = sqrt(t1)
            If (t1 .Gt. forcemax) forcemax = t1
         End Do
      End Do

      Deallocate (rfmt, grfmt)
      Call timesec (ts1)
      timefor = timefor + ts1 - ts0
!      write(*,*) 'force'
      Return
End Subroutine
!EOC
