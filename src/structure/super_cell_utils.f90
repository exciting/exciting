
module super_cell_utils
   use precision, only: dp
   use asserts, only: assert

   implicit none
   private
   public :: get_translation_vectors, &
             supercell_atomic_positions, &
             TranslationIntegers_type, &
             extend_rmt_natoms

   !> Translation integers type containing methods to get
   !> total number of translations and checking if their
   !> values are valid.
   type TranslationIntegers_type
      integer :: i(2)
      integer :: j(2)
      integer :: k(2)
   contains
      !> Get the total number of integers (hence translations)
      procedure :: total => get_total_number_of_integers
      !> Check if integer values are valid loop limits. 
      procedure :: integers_valid => check_integers_valid
   end type

contains

   !> Calculates total number of translations.
   function get_total_number_of_integers(this) result(total)
      class(TranslationIntegers_type), intent(inout) :: this
      integer :: total

      total = (this%i(2) - this%i(1) + 1)* &
              (this%j(2) - this%j(1) + 1)* &
              (this%k(2) - this%k(1) + 1)

   end function

   !> Checks if class arguments are valid loop limits.
   function check_integers_valid(this) result(integers_valid)
      class(TranslationIntegers_type), intent(inout) :: this
      logical :: integers_valid

      integers_valid = this%i(2) >= this%i(1) .and. &
                       this%j(2) >= this%j(1) .and. &
                       this%k(2) >= this%k(1)

   end function

   !> Generates a set of translation vectors.
   !>
   !> Given lattice vectors \(\mathbf{L}\) in columnwise format, and a set of integers \(\mathbf{n}\),
   !> generates a set of translation vectors \(\mathbf{T}\), where \(\mathbf{T}\) is defined as:
   !>
   !>  \[ \mathbf{T} = \mathbf{L}\mathbf{n} \]
   function get_translation_vectors(lattice_vect, n) result(translation)

      !> Lattice vectors
      real(dp), intent(in) :: lattice_vect(3, 3)
      !> Translation vectors
      type(TranslationIntegers_type), intent(inout) :: n

      integer :: itrans, i1, i2, i3

      !> Translated vector
      real(dp), allocatable:: translation(:, :)
      allocate(translation(3, n%total()))

      call assert(n%integers_valid(), message='For each dimension the second integer &
                  value has to be larger than or equal to the first integer value.')

      itrans = 0
      do i1 = n%i(1), n%i(2)
         do i2 = n%j(1), n%j(2)
            do i3 = n%k(1), n%k(2)
               itrans = itrans + 1
               translation(:, itrans) = matmul(lattice_vect, [i1, i2, i3])
            end do
         end do
      end do

   end function get_translation_vectors

   !> Given a set of translation vectors \( \mathbf{T} \), containing the vectors pointing 
   !> to all cells, and the basis \(\mathbf{atpos}\), will generate all atomic positions 
   !> \( \mathbf{r} \) for each unit cell.
   !> \[ 
   !>    \mathbf{r} = \mathbf{T} + \mathbf{atposc}
   !> \]
   function supercell_atomic_positions(translation, atomic_postion, natoms) result(positions)

      !> Translation vectors
      real(dp), intent(in) :: translation(:, :)
      !> Atomic positions
      real(dp), intent(in) :: atomic_postion(:, :, :)
      !> Number of atoms per species
      integer, intent(in) :: natoms(:)

      integer :: iatom, i, is, ia
      real(dp), allocatable :: t_vector(:)

      !> Supercell atomic positions
      real(dp), allocatable :: positions(:, :)
      allocate (positions(3, sum(natoms)*size(translation, 2)))

      call assert(size(translation, 1) == 3, message= &
                  "Expected vector component along the first dimension is three.")
      call assert(size(atomic_postion, 1) == 3, message= &
                  "Expected vector component along the first dimension is three.")

      iatom = 0
      do i = 1, size(translation, 2)
         t_vector = translation(:, i)
         do is = 1, size(natoms)
            do ia = 1, natoms(is)
               iatom = iatom + 1
               positions(:, iatom) = t_vector + atomic_postion(:, ia, is)
            end do
         end do
      end do

   end function supercell_atomic_positions

  
   !> Given the muffin-tin radius \(\mathbf{rmt}(is)\), for each species \(is\), and the number of atoms per 
   !> species \(\mathbf{natoms}(is)\), expands the \(\mathbf{rmt}\) vector of size \(ns\), the total number
   !> of species, to the size of the total number of atoms in one cell \(na\), containing the corresponding 
   !> rmt value for each atom instead of for each species.
   function extend_rmt_natoms(natoms, rmt) result(extended_rmt)
      !> Number of atoms per species
      integer, intent(in) :: natoms(:)
      !> Muffin-tin radius for each species
      real(dp), intent(in) :: rmt(:)

      integer :: is, ia, iatom, ns

      !> Muffin-tin radius for each atom
      real(dp), allocatable :: extended_rmt(:)

      call assert(size(natoms) == size(rmt), message= &
                  "Dimensions of natoms and rmt have to be equal.") 
      allocate(extended_rmt(sum(natoms)))

      ns = size(natoms)
      iatom = 0
      do is = 1, ns
         do ia = 1, natoms(is)
            iatom = iatom + 1
            extended_rmt(iatom) = rmt(is)
         end do
      end do

   end function extend_rmt_natoms


end module
