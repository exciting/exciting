!> @file exciting_idiel_interface.F90
!> @brief Safe interface to the IDieL library within the exciting code.
!>
!> This module offers a wrapper around the IDieL library for computing
!> anisotropic averages, while providing fallback behavior
!> when IDieL is not compiled in.

module exciting_idiel_interface

#if defined(IDIEL)
    use idiel, only: idiel_t
#endif

    use precision, only: i32, dp
    use modinput, only: scrcoul_type

    implicit none

    private

    public :: init_idiel_handler, destroy_idiel_handler, set_head_wings_idiel, get_head_wings_idiel, invert_body, &
              idiel_2d_anisotropic_average_idiel, idiel_3d_anisotropic_average_idiel, get_body_idiel, free_memory_idiel

#if defined(IDIEL)
    !> The IDieL library handler
    type(idiel_t), target :: idiel_handler
#endif

    type(scrcoul_type), pointer :: scrcoul
    logical :: hermitian
    logical :: idiel_averaging

contains

   !> Initializes the IDieL handler and required data structures.
   !>
   !> Sets up system geometry, atomic species, and reduced positions,
   !> and configures the dielectric averaging mode based on input.
   !>
   !> If compiled without IDieL and averaging mode is 'anisotropic'
   !> or 'anisotropic-2d' terminates execution with an error.
   subroutine init_idiel_handler()
      use modinput,     only: scrcoul_type, input
      use mod_atoms,    only: nspecies, natoms, idxas, natmtot
      use mod_lattice,  only: avec
      use modinput,     only: input
      use mod_device_offload,    only: device_world
      use modmpi, only: terminate
      use modgw,  only: freq

      real(dp), allocatable      :: reduced_pos(:,:)
      integer(i32), allocatable  :: elements(:)
      integer(i32) :: iatom, is, ia, system_dimension

      scrcoul    => input%gw%scrcoul
      idiel_averaging  = (trim(scrcoul%averaging) == 'anisotropic') .or. &
                   (trim(scrcoul%averaging) == 'anisotropic-2d')

      !> For purely imaginary frequencies the RPA dielectric matrix is Hermitian
      !> so the upper wing is the complex conjugate of the lower wing
      !> which permits a more efficient algorithm 
      !> for the averaging
      hermitian = (freq%fconv == 'imfreq')

#if defined(IDIEL)
      allocate(reduced_pos(natmtot, 3))
      allocate(elements(natmtot))

      ! Get the system dimensionality
      system_dimension = merge(2, 3, trim(scrcoul%averaging) == 'anisotropic-2d')

      ! We need to write the reduced atomic positions in an IDieL 
      ! friendly format
      iatom = 1
      do is = 1, nspecies
         do ia = 1, natoms(is)
            reduced_pos(iatom,:) =  &
               input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
            ! In exciting atomic numbers are negative 
            elements(iatom) = iabs(input%structure%speciesarray(is)%species%atomicNumber)
            iatom = iatom + 1
         end do
      end do
      ! Call the handler creation
      ! Note that in exciting the lattice vectors are given 
      ! in columns, and IDieL wants them in rows.
      call idiel_handler%init_common(transpose(avec), transpose(reduced_pos), &
            elements, input%gw%ngridq, dim=system_dimension, device_world=device_world)
#else
      if ( idiel_averaging ) then
          call terminate("Error scrcoul%averaging  anisotropic or anisotropic-2d require IDieL support")
      end if
#endif
   end subroutine init_idiel_handler

   !> Frees all resources used by the IDieL handler.
   subroutine destroy_idiel_handler()
#if defined(IDIEL)
      call idiel_handler%clean()
#endif
   end subroutine destroy_idiel_handler

   !> Sets the head and wings of the dielectric matrix in IDieL, for a given frequency
   !> at the Gamma point
   subroutine set_head_wings_idiel(head, wingU, wingL, freq_idx, is_Gamma)
       !> Head of the dielectric matrix, the last index is the frequency index
       complex(dp), intent(in)  :: head(:,:,:)
       !> Upper wing of the dielectric matrix, the last index is the frequency index
       complex(dp), intent(in)  :: wingU(:,:,:)
       !> Lower wing of the dielectric matrix, the last index is the frequency index
       complex(dp), intent(in)  :: wingL(:,:,:)
       !> Frequency index
       integer(i32), intent(in) :: freq_idx
       !> Are we dealing with the Gamma point?
       logical,      intent(in) :: is_Gamma
#if defined(IDIEL)
       if (is_Gamma .and. idiel_averaging) call idiel_handler%set_dielectric_blocks(head(:,:,freq_idx), wingU(:,:,freq_idx), wingL(:,:,freq_idx))
#endif
   end subroutine set_head_wings_idiel

   !> Retrieves the head and wings of the dielectric matrix from IDieL, for a given frequency 
   !> at Gamma
   subroutine get_head_wings_idiel(head, wingU, wingL, freq_idx, is_Gamma)
       !> Head of the dielectric matrix, the last index is the frequency index
       complex(dp), intent(inout)  :: head(:,:,:)
       !> Upper wing of the dielectric matrix, the last index is the frequency index
       complex(dp), intent(inout)  :: wingU(:,:,:)
       !> Lower wing of the dielectric matrix, the last index is the frequency index
       complex(dp), intent(inout)  :: wingL(:,:,:)
       !> Frequency index
       integer(i32), intent(in)    :: freq_idx
       !> Are we dealing with the Gamma point?
       logical,      intent(in)    :: is_Gamma
#if defined(IDIEL)
       if (is_Gamma .and. idiel_averaging) then
           head(1,1,freq_idx)   = idiel_handler%idiel_head
           ! The wings are zero because they are odd functions of q,
           ! and thus their average vanishes when integrated symmetrically around q=0.
           wingL(:,1,freq_idx)  = idiel_handler%idiel_wingL
           wingU(:,1,freq_idx)  = idiel_handler%idiel_wingU
       end if
#endif
   end subroutine get_head_wings_idiel

   !> Inverts the body of the dielectric matrix.
   subroutine invert_body(body)
       use inverse, only: invert_LU
       !> The body to compute the inverse of
       complex(dp), pointer, intent(inout)  :: body(:,:)
#if defined(IDIEL)
        !> Invert the body using IDieL; as it 
        !> allows for GPU offload of the inversion
        call idiel_handler%invert_body(body)
        body => idiel_handler%Binv_data
#else
        call invert_LU(body)
#endif
   end subroutine invert_body

   !> Retrieves the inverted dielectric matrix body block from IDieL at Gamma.
   subroutine get_body_idiel(body, is_Gamma)
       !> The body of the dielectric matrix
       complex(dp), pointer, intent(inout)  :: body(:,:)
       !> Are we in the Gamma point
       logical, intent(in) :: is_Gamma
#if defined(IDIEL)
        if (is_Gamma .and. idiel_averaging) then
            body => idiel_handler%idiel_body
        end if
#endif
   end subroutine get_body_idiel

   !> @brief Computes anisotropic average of the dielectric matrix for 2D systems.
   subroutine idiel_2d_anisotropic_average_idiel()
        use scrcoul_low_dim, only: set_singc12
        use constants, only: zone
#if defined(IDIEL)
        integer(i32) :: i, nb

        ! For cases where the RPA dielectric matrix is Hermitian,
        ! more efficient algorithms are used.
        ! See the assignment of the 'hermitian' variable for more information.
        call idiel_handler%compute_anisotropic_avg_scrcoulomb_2d(hermitian)
        call set_singc12
        idiel_handler%idiel_head = idiel_handler%idiel_head  + zone
        nb = size(idiel_handler%idiel_body,1)
        do i = 1, nb
            idiel_handler%idiel_body(i,:) = idiel_handler%idiel_body(i,:) + zone
        end do
#endif
   end subroutine idiel_2d_anisotropic_average_idiel

   !> Computes anisotropic average of the dielectric matrix for 3D systems.
   subroutine idiel_3d_anisotropic_average_idiel()
#if defined(IDIEL)
        ! For cases where the RPA dielectric matrix is Hermitian,
        ! more efficient algorithms are used.
        ! See the assignment of the 'hermitian' variable for more information.
        call idiel_handler%compute_anisotropic_avg_inversedielectric_3d(hermitian)
#endif
   end subroutine idiel_3d_anisotropic_average_idiel

   !> Free the body memory inside IDieL
   subroutine free_memory_idiel()
#if defined(IDIEL)
      if (allocated(idiel_handler%Binv_data))   deallocate(idiel_handler%Binv_data)
      if (allocated(idiel_handler%idiel_body))  deallocate(idiel_handler%idiel_body)
      if (allocated(idiel_handler%idiel_wingL)) deallocate(idiel_handler%idiel_wingL)
      if (allocated(idiel_handler%idiel_wingU)) deallocate(idiel_handler%idiel_wingU)
#endif
   end subroutine free_memory_idiel 

end module exciting_idiel_interface

