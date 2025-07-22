module invert_dielectric_function
   use asserts, only: assert
   use precision, only: dp, i32
   use constants, only: zone
   use exciting_idiel_interface, only: set_head_wings_idiel, get_head_wings_idiel, invert_body, &
                                       idiel_2d_anisotropic_average_idiel, free_memory_idiel, &
                                       idiel_3d_anisotropic_average_idiel, get_body_idiel
   implicit none
   private
   public :: calcinveps

contains
   !> Compute the inverse of the dielectric function, including
   !> applying spherical averaging.
   subroutine calcinveps(iomstart, iomend, gamma, scrcoul, freqtype, symt2, epsilon, epsw1, epsw2, epsh, eps00, time_dfinv)
      use modinput, only: scrcoul_type
      use inverse, only: invert_LU
      !>Starting frequency index
      integer(i32), intent(in) :: iomstart
      !>Ending frequency index
      integer(i32), intent(in) :: iomend
      !>Is true if gamma point is present
      logical, intent(in) :: gamma
      !>Input scrcoul element 
      type(scrcoul_type), optional, intent(in) :: scrcoul
      !>Select whether frequencies are real or imaginary
      character(6), intent(in) :: freqtype
      !>Symmetrization tensor
      real(dp), optional, intent(in) :: symt2(3, 3, 3, 3)
      !>In: Body of the dielectric function. Out: Body of the inverse dielectric function
      complex(dp), target, intent(inout) :: epsilon(:, :, :)
      !>In: 1st wing of the dielectric function. Out: 1st wing of the inverse dielectric function
      complex(dp), optional, intent(inout) :: epsw1(:, :, :)
      !>In: 2nd wing of the dielectric function. Out: 2nd wing of the inverse dielectric function
      complex(dp), optional, intent(inout) :: epsw2(:, :, :)
      !>In: Head of the dielectric function. Out: Head of the inverse dielectric function
      complex(dp), optional, intent(inout) :: epsh(:, :, :)
      !>Symmetrized dielectric tensor
      complex(dp), optional, intent(out) :: eps00(:, :, :)
      !>Timing
      real(dp), intent(inout) :: time_dfinv


      integer(i32) :: iom, mbsiz, i, j
      integer(i32) :: im
      integer(i32) :: info, lwork
      real(dp)    :: tstart, tend
      complex(dp), pointer, contiguous :: eps(:,:)

      character(len=10), parameter :: sname="calcinveps"

      call timesec(tstart)
      if( gamma ) then
        call assert( present(scrcoul) .and. present(symt2) .and. present(epsw1) .and. &
          present(epsw2) .and. present(epsh) .and. present(eps00), &
          'Optional arguments must be present when gamma is true' )
      end if
      mbsiz = size(epsilon, 1)

      ! lopp over frequencies
      do iom = iomstart, iomend

         ! array for body and its inverse
         eps => epsilon(1:mbsiz,1:mbsiz,iom)

         ! We invert the body; if IDieL is available
         ! we use it as it allows device offload
         ! If using IDieL this reassociates the pointer
         ! to the interal buffer holding of IDieL holding
         ! the actual inverse of the body.
         call invert_body(eps)

         !averaging of eps for q->0
         if (Gamma) then
            call set_head_wings_idiel(epsh, epsw1, epsw2, iom, Gamma)
            call angular_averaging(iom, symt2, scrcoul, eps, epsw1, epsw2, epsh, eps00)
            call get_head_wings_idiel(epsh, epsw1, epsw2, iom, Gamma)
            call get_body_idiel(eps, Gamma)
            epsh(1,1,iom) = epsh(1,1,iom) - zone !\epsilon^{-1}_{00}-1
         endif

         ! If using IDieL `eps` has been reassociated to
         ! internal buffers of IDieL, thus is not associated
         ! anymore with the global `epsilon`. Therefore
         ! we need to update it back.
         ! Some compilers create for some reason a temporary to make the 
         ! copy if one writes epsilon(1:mbsiz,1:mbsiz,iom) = eps(1:mbsiz,1:mbsiz)
         ! The same does not happen with the explicit do copy.
         !$omp parallel do collapse(2) default(none) private(i,j) shared(epsilon,eps,mbsiz,iom)
         do i = 1, mbsiz
           do j = 1, mbsiz
             epsilon(j,i,iom) = eps(j,i)
           end do
         end do
         !$omp end parallel do

         ! Update diagonal: epsilon^{-1}_{ij} - \delta_{ij}
         do im = 1, mbsiz
            epsilon(im,im,iom) = epsilon(im,im,iom)-zone
         end do

         ! Nullify the pointer
         nullify(eps)

         ! If using IDieL free memory
         call free_memory_idiel()

      end do

      call timesec(tend)
      time_dfinv = time_dfinv+tend-tstart
   end subroutine

   !>The averaging schemes mainly follow:
   !>C. Freysoldt, P. Eggert, P. Rinke, A. Schindlmayr, R. W. Godby, and M. Scheffler, Comput. Phys. Commun. 176, 1 (2007)
   !>In the limit of \( q\to 0\) the inverse dielectric function can be written as:
   !>\begin{align}
   !>\epsilon^{-1}_{00} &= \frac{1}{\hat{\mathbf{q}}L\hat{\mathbf{q}}} \\
   !>\epsilon^{-1}_{\mu 0} &= -\epsilon^{-1}_{00}\mathbf{s}_\mu\cdot\hat{\mathbf{q}} \\
   !>\epsilon^{-1}_{0 \mu} &= -\epsilon^{-1}_{00}\mathbf{t}_\mu\cdot\hat{\mathbf{q}} \\
   !>\epsilon^{-1}_{\mu\nu} &= B^{-1}_{\mu\nu} + \epsilon^{-1}_{00}(\mathbf{s}_\mu\cdot\hat{\mathbf{q}})(\mathbf{t}_\nu\cdot\hat{\mathbf{q}}) \\
   !>\end{align}
   !>where \( \hat{\mathbf{q}}\) is the direction in which the limit is taken, \( L \) a 3x3 tensor and \( \mathbf{s}_\mu, \mathbf{t}_\mu \) are vectors.
   subroutine angular_averaging(iom, symt2, scrcoul, eps, epsw1, epsw2, epsh, eps00)
      use modinput, only: scrcoul_type
      use scrcoul_low_dim, only: apply_2d_limit
      use mod_lattice, only: bvec
      use mod_coulomb_potential, only: rcut
      use modmpi, only: terminate
      use ieee_arithmetic, only: ieee_value, ieee_signaling_nan
      !>Frequency index
      integer(i32), intent(in) :: iom
      !>Symmetrization tensor
      real(dp), intent(in) :: symt2(3, 3, 3, 3)
      !>Input scrcoul element 
      type(scrcoul_type), intent(in) :: scrcoul
      !>In: Inverse of the body of the dielectric function. Out: Body of the inverse dielectric function
      complex(dp), intent(inout) :: eps(:, :)
      !>In: 1st wing of the dielectric function. Out: 1st wing of the inverse dielectric function
      complex(dp), intent(inout) :: epsw1(:, :, :)
      !>In: 2nd wing of the dielectric function. Out: 2nd wing of the inverse dielectric function
      complex(dp), intent(inout) :: epsw2(:, :, :)
      !>In: Head of the dielectric function. Out: Head of the inverse dielectric function
      complex(dp), intent(inout) :: epsh(:, :, :)
      !>Symmetrized dielectric tensor
      complex(dp), intent(out) :: eps00(:, :, :)

      select case(trim(scrcoul%averaging))
         case("isotropic")
            call isotropic_averaging(iom, symt2, scrcoul%q0eps, eps, epsw1, epsw2, epsh, eps00)
         case("2d")
            call apply_2d_limit(bvec, rcut, symt2, eps(:, :), epsh(:, :, iom), epsw1(:, :, iom), epsw2(:, :, iom))
            !Warning Hack:
            epsh(1,1,iom) = epsh(1,1,iom)+zone
         case("anisotropic-2d")
            call idiel_2d_anisotropic_average_idiel()
         case("anisotropic")
            call idiel_3d_anisotropic_average_idiel()
         case default
            call terminate("FATAL ERROR(angular_averaging): scrcoul%averaging is in an invalid state. " // &
                           "This condition should be unreachable under correct program logic. " // &
                           "Reaching this point indicates a critical internal failure.")
      end select

      ! After averaging, the head is no longer a tensor but a scalar. However, the global array is still reused, with only the (1,1)
      ! element being updated - the rest of the elements keep their previous values.
      ! For task-based GW, where all elements of the head are printed in the inverse epsilon file.
      ! This can be confusing, since only the (1,1) element is actually meaningful after averaging.
      ! To make this clear, we set the unused elements to signaling NaNs.
      ! That way, if anyone tries to use them in arithmetic operations, it will raise an error;
      ! making it obvious that those parts of the array are now garbage and should not be used.
      ! The same logic applies to the wings.
      associate ( complex_nan => cmplx(ieee_value(0.0_dp, ieee_signaling_nan), &
                                       ieee_value(0.0_dp, ieee_signaling_nan), kind=dp))

          epsh(2:3, 2:3, iom) = complex_nan
          epsh(2:3, 1,   iom) = complex_nan
          epsh(1,   2:3, iom) = complex_nan
          epsw1(:, 2:3, iom)  = complex_nan
          epsw2(:, 2:3, iom)  = complex_nan

      end associate

   end subroutine


   !>Calculates the symmetrised dielectric tensor, and inverse of wings 1 and 2.
   subroutine symmetrised_dielectric_tensor(iom, eps, epsw1, epsw2, epsh, symt2, L, s, t)
      use xlapack, only: matrix_multiply
      !>Frequency index
      integer, intent(in) :: iom
      !>In: Inverse of the body of the dielectric function.
      complex(dp), intent(in) :: eps(:, :)
      !>In: 1st wing of the dielectric function.
      complex(dp), intent(in) :: epsw1(:, :, :)
      !>In: 2nd wing of the dielectric function.
      complex(dp), intent(in) :: epsw2(:, :, :)
      !>In: Head of the dielectric function. 
      complex(dp), intent(in) :: epsh(:, :, :)
      !>Symmetrization tensor
      real(dp), intent(in) :: symt2(3, 3, 3, 3)
      !>Symmetrized dielectric tensor
      complex(dp), intent(out) :: L(3, 3)
      !>Vector that defines the 1st wing of the inverse dielectric function
      complex(dp), allocatable, intent(out) :: s(:, :)
      !>Vector that defines the 2nd wing of the inverse dielectric function
      complex(dp), allocatable, intent(out) :: t(:, :)
      
      integer :: mbsiz, iop, jop
      complex(dp) :: dtns(3,3), lfe(3,3)
      complex(dp), allocatable :: u(:,:), v(:,:)

      mbsiz = size(eps, 1)
      ! definition of U_{\alpha} (column wing) (B.7)
      allocate(u(mbsiz,3))
      u(:,:) = epsw1(:,:,iom)
      ! row wing
      allocate(v(3,mbsiz))
      v(:,:) = transpose(epsw2(:,:,iom))

      ! definition of S{\alpha} (B.13)
      allocate(s(mbsiz,3))
      call matrix_multiply(eps, u, s)
      allocate(t(3,mbsiz))
      call matrix_multiply(v, eps, t)

     
      ! definition of L = H-conjg(U)*S (B.14) = \epsilon_{00}+LFE
      call matrix_multiply(v, s, lfe)
      L(:,:) = epsh(:,:,iom) - lfe

      ! symmetrize the dielectric tensor

      dtns = L
      do iop = 1, 3
         do jop = 1, 3
            call symt2app(iop, jop, 1, symt2, dtns, L(iop,jop))
         end do
      end do

      deallocate(u, v)

   end subroutine

   !>In this routine one specific direction for \( \hat{\mathbf{q}} \) is passed as the input parameter `q0eps`.
   subroutine isotropic_averaging(iom, symt2, q0eps, eps, epsw1, epsw2, epsh, eps00)
      use asserts, only: assert
      use math_utils, only: all_zero
      !>Frequency index
      integer(i32), intent(in) :: iom
      !>Symmetrization tensor
      real(dp), intent(in) :: symt2(3, 3, 3, 3)
      !>Direction in which the limit is taken
      real(dp), intent(in) :: q0eps(3)
      !>In: Inverse of the body of the dielectric function. Out: Body of the inverse dielectric function
      complex(dp), intent(inout) :: eps(:, :)
      !>In: 1st wing of the dielectric function. Out: 1st wing of the inverse dielectric function
      complex(dp), intent(inout) :: epsw1(:, :, :)
      !>In: 2nd wing of the dielectric function. Out: 2nd wing of the inverse dielectric function
      complex(dp), intent(inout) :: epsw2(:, :, :)
      !>In: Head of the dielectric function. Out: Head of the inverse dielectric function
      complex(dp), intent(inout) :: epsh(:, :, :)
      !>Symmetrized dielectric tensor
      complex(dp), intent(out) :: eps00(:, :, :)

      ! local variables
      integer(i32) :: mbsiz
      integer(i32) :: j1, j2, i
      integer(i32) :: iop, jop
      real(dp) :: q0eps_dot_q0eps
      complex(dp) :: L(3,3), L_diag(3), dtns(3,3)
      complex(dp), allocatable :: s(:,:), t(:,:)
      !> Tolerance for zero
      real(dp), parameter :: tol = 1.e-8_dp

      call assert(.not. all_zero(q0eps), "q0eps should not be zero")

      mbsiz = size(eps, 1)
      call symmetrised_dielectric_tensor(iom, eps, epsw1, epsw2, epsh, symt2, L, s, t)

      !-------------------------------------------
      ! Store the symmetrized macroscopic tensor
      !-------------------------------------------
      eps00(:,:,iom) = L(:,:)

      L_diag = [(L(i, i), i = 1, 3)]

      call assert(maxval(abs(L_diag))>tol, "Diagonal elements of dielectric tensor should not be zero")
      
      !====================================================
      ! calculate the averaged inverse dielectric function
      !====================================================

      q0eps_dot_q0eps = dot_product(q0eps, q0eps)
         
      epsh(1, 1, iom) = q0eps_dot_q0eps / dot_product(L_diag, q0eps)

      epsw1(:,1,iom) = -epsh(1,1,iom)* &
      &                (s(:,1)*q0eps(1)+s(:,2)*q0eps(2)+s(:,3)*q0eps(3))/dsqrt(q0eps_dot_q0eps)

      epsw2(:,1,iom) = -epsh(1,1,iom)* &
      &                (t(1,:)*q0eps(1)+t(2,:)*q0eps(2)+t(3,:)*q0eps(3))/dsqrt(q0eps_dot_q0eps)

      do j2 = 1, mbsiz
         do j1 = 1, mbsiz
            eps(j1,j2) = eps(j1,j2) + &
            &           epsw1(j1,1,iom)*epsw2(j2,1,iom)/epsh(1,1,iom)
         end do
      end do

      deallocate(s, t)
   end subroutine

end module
