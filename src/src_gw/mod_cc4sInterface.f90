!> This module contains procedures to store and manipulate properties needed as input
!> for the coupled cluster code cc4s
module mod_cc4sInterface
    use constants, only: zzero
    use gw_io, only: build_file_name, write_to_file, read_from_file
    use modmpi, only: terminate
    use precision, only: i32, dp
#include "offload.fpp"

    implicit none

    private

    !> Coulomb Vertex
    complex(dp), protected, allocatable :: coulomb_vertex(:,:,:)

    !> cc4s input files
    character(len=*), parameter :: file_name_coulomb_vertex = 'CoulombVertex.elements'
    character(len=*), parameter :: yaml_name_coulomb_vertex = 'CoulombVertex.yaml'
    character(len=*), parameter :: file_name_eigen_energies = 'EigenEnergies.elements'
    character(len=*), parameter :: yaml_name_eigen_energies = 'EigenEnergies.yaml'
    character(len=*), parameter :: yaml_name_orbital_properties = 'State.yaml'


    public :: init_coulomb_vertex, delete_coulomb_vertex, &
             compute_coulomb_vertex, prepare_dft_eigenvalues, &
             write_coulomb_vertex_info_to_yaml, write_coulomb_vertex_to_file, &
             write_eigenenergies_to_yaml, write_scf_energies_to_file, &
             write_orbital_properties_to_yaml, &
             file_name_coulomb_vertex, yaml_name_coulomb_vertex, &
             file_name_eigen_energies, yaml_name_eigen_energies, &
             yaml_name_orbital_properties, coulomb_vertex, &
             write_test_output

    
contains

subroutine init_coulomb_vertex(mbsiz, nstates)
  use constants, only: zzero
  implicit none
  !> Mixed basis size
  integer(i32), intent(in) :: mbsiz
  !> Number of states
  integer(i32), intent(in) :: nstates

  ! Clean previous
  OMP_OFFLOAD target exit data map(always, delete: coulomb_vertex) if(allocated(coulomb_vertex))
  if (allocated(coulomb_vertex)) deallocate(coulomb_vertex)

  ! allocate coulomb_vertex
  allocate(coulomb_vertex(mbsiz+1, nstates, nstates), source=zzero)
  OMP_OFFLOAD target enter data map(always, to: coulomb_vertex)

end subroutine init_coulomb_vertex


subroutine delete_coulomb_vertex()

  if (allocated(coulomb_vertex)) then
    OMP_OFFLOAD target exit data map(always, delete: coulomb_vertex)
    deallocate(coulomb_vertex)
  end if

end subroutine delete_coulomb_vertex


subroutine write_coulomb_vertex_info_to_yaml(dim_2_3, dim_1, e_unit, file_name)
  implicit none

  !> Size of dim 2 and 3 of Coulom vertex (Number of states)
  integer(i32), intent(in) :: dim_2_3
  !> Size of dim 1 of Coulom vertex (mbsiz+1)
  integer(i32), intent(in) :: dim_1
  !> Energy unit used
  real(dp),  intent(in) :: e_unit
  !> Name of the output file
  character(len=*), intent(in) :: file_name

  !local
  integer(i32) :: file_unit

  open (newunit=file_unit, file=file_name)
  write (file_unit, '(A12)') "version: 100"
  write (file_unit, '(A12)') "type: Tensor"
  write (file_unit, '(A21)') "scalarType: Complex64"
  write (file_unit, '(A11)') "dimensions:"
  write (file_unit, '(A9, 1X, I0)') "- length:", dim_1
  write (file_unit, '(2X, A20)') "type: AuxiliaryField"
  write (file_unit, '(A9, 1X, I0)') "- length:", dim_2_3
  write (file_unit, '(2X, A11)') "type: State"
  write (file_unit, '(A9, 1X, I0)') "- length:", dim_2_3
  write (file_unit, '(2X, A11)') "type: State"

  write (file_unit, '(A9)') "elements:"
  write (file_unit, '(2X, A20)') "type: IeeeBinaryFile"
  write (file_unit, '(A5, 1X, F12.7)') "unit:", e_unit

  write (file_unit, '(A9)') "metaData:"
  write (file_unit, '(2X, A11)') "halfGrid: 0"

  close (unit=file_unit)

end subroutine write_coulomb_vertex_info_to_yaml


!> Write the (q-dependent) Coulomb vertex to a binary file.
!>
!> When compiled with MPI, the last dimension of coulomb_vertex is
!> distributed in contiguous blocks across `nproc` ranks (any remainder
!> elements go to the last rank), and each rank performs a SINGLE
!> collective write into its slice of the file, located via an
!> MPI subarray file view. When compiled without MPI, the whole array
!> is written directly with one serial stream write.
!>
!> NOTE: `iq` is appended to `file_name` so repeated calls for
!> different q-points don't overwrite each other. cc4s can not
!> handle multiple q-points yet, so you need to choose the 
!> file corresponding to the q-point you would like to use
!> and rename it to CoulombVertex.elements .
subroutine write_coulomb_vertex_to_file(iq, coulomb_vertex, file_name)

  use iso_fortran_env, only: int64

#ifdef MPI
  use mpi_f08, only: MPI_Comm, MPI_Datatype, MPI_File, MPI_OFFSET_KIND, &
       MPI_DOUBLE_COMPLEX, MPI_INFO_NULL, MPI_MODE_CREATE, MPI_MODE_WRONLY, &
       MPI_SUCCESS, MPI_COMM_WORLD, MPI_ORDER_FORTRAN, MPI_STATUS_IGNORE, &
       MPI_Comm_size, MPI_Comm_rank, MPI_Cart_create, MPI_Cart_coords, &
       MPI_Comm_free, MPI_Type_create_subarray, MPI_Type_contiguous, &
       MPI_Type_commit, MPI_Type_free, &
       MPI_File_open, MPI_File_set_view, MPI_File_write_at_all, MPI_File_close, &
       MPI_Abort
#endif

  implicit none

  !> current q-point index (right now, always 1)
  integer(i32), intent(in) :: iq
  !> Coulomb vertex
  complex(dp), intent(in)  :: coulomb_vertex(:, :, :)
  !> Name of the output file
  character(len=*), intent(in) :: file_name

  character(len=300) :: fname_buffer
  character(len=:), allocatable :: full_filename
  integer :: ierr

#ifdef MPI
  integer, parameter :: ndims = 3
  integer :: dims(ndims), coords(ndims)
  integer :: int_shape(ndims), int_subshape(ndims), int_substart(ndims)
  logical :: periods(ndims) = .false.
  logical :: reorder = .false.
  type(MPI_Comm)     :: comm
  type(MPI_Datatype) :: subarr, plane_type
  type(MPI_File)     :: fh
  integer :: myrank, nproc, blocklen, planecount
  integer(MPI_OFFSET_KIND) :: disp
  integer(int64) :: i1, i2, i3
  integer(int64) :: vertex_shape(ndims), subvertex_shape(ndims), subvertex_start_coord(ndims)
  complex(dp), allocatable :: subvertex(:, :, :)
#else
  integer :: funit
#endif

  write(fname_buffer, '(A,A,I0)') trim(file_name), '.q', iq
  full_filename = trim(fname_buffer)

#ifdef MPI

  i1 = size(coulomb_vertex, 1, kind=int64)
  i2 = size(coulomb_vertex, 2, kind=int64)
  i3 = size(coulomb_vertex, 3, kind=int64)
  vertex_shape = [i1, i2, i3]

  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  ! 1-D Cartesian communicator: parallelise only over the last
  ! dimension of coulomb_vertex.
  dims = [1, 1, nproc]
  call MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, comm, ierr)
  call MPI_Comm_rank(comm, myrank, ierr)
  call MPI_Cart_coords(comm, myrank, ndims, coords, ierr)

  ! Distribute i3 as evenly as possible; remainder goes to last rank.
  subvertex_shape    = vertex_shape
  subvertex_shape(3) = i3 / nproc
  subvertex_start_coord    = 0_int64
  subvertex_start_coord(3) = coords(3) * (i3 / nproc)
  if (myrank == nproc - 1) then
     subvertex_shape(3) = subvertex_shape(3) + mod(i3, int(nproc, int64))
  end if

  allocate(subvertex(i1, i2, subvertex_shape(3)))
  subvertex = coulomb_vertex(:, :, subvertex_start_coord(3) + 1 : &
                                    subvertex_start_coord(3) + subvertex_shape(3))

  ! MPI_Type_create_subarray takes default-kind integer arrays.
  int_shape    = int(vertex_shape,          kind=kind(1))
  int_subshape = int(subvertex_shape,       kind=kind(1))
  int_substart = int(subvertex_start_coord, kind=kind(1))

  call MPI_Type_create_subarray(ndims, int_shape, int_subshape, int_substart, &
       MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, subarr, ierr)
  call MPI_Type_commit(subarr, ierr)

  disp = 0_MPI_OFFSET_KIND
  call MPI_File_open(comm, full_filename, MPI_MODE_WRONLY + MPI_MODE_CREATE, &
       MPI_INFO_NULL, fh, ierr)
  if (ierr /= MPI_SUCCESS) then
     write(*,*) "write_coulomb_vertex_to_file: could not open ", trim(full_filename)
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if

  call MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, subarr, 'native', MPI_INFO_NULL, ierr)

  if (i1 * i2 > huge(blocklen) .or. subvertex_shape(3) > huge(planecount)) then
     write(*,*) "write_coulomb_vertex_to_file: plane size or plane count", &
                " exceeds default INTEGER range on rank", myrank
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if
  blocklen   = int(i1 * i2, kind=kind(1))
  planecount = int(subvertex_shape(3), kind=kind(1))

  call MPI_Type_contiguous(blocklen, MPI_DOUBLE_COMPLEX, plane_type, ierr)
  call MPI_Type_commit(plane_type, ierr)

  call MPI_File_write_at_all(fh, 0_MPI_OFFSET_KIND, subvertex, planecount, &
       plane_type, MPI_STATUS_IGNORE, ierr)
  if (ierr /= MPI_SUCCESS) then
     write(*,*) "write_coulomb_vertex_to_file: collective write failed on rank", myrank
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if

  call MPI_Type_free(plane_type, ierr)

  call MPI_File_close(fh, ierr)
  call MPI_Type_free(subarr, ierr)
  call MPI_Comm_free(comm, ierr)
  deallocate(subvertex)

#else

  ! Serial fallback: no MPI needed, just stream the whole array out.
  open(newunit=funit, file=full_filename, form='unformatted', access='stream', &
       status='replace', action='write', iostat=ierr)
  if (ierr /= 0) then
     call terminate('write_coulomb_vertex_to_file: could not open ' // trim(full_filename))
  end if
  write(funit) coulomb_vertex
  close(funit)

#endif

end subroutine write_coulomb_vertex_to_file

!> Writes information on and the eigenenergies itself to an output file. 
subroutine write_eigenenergies_to_yaml(e_unit, e_fermi, scf_energies, file_name)
  !> Energy unit
  real(dp), intent(in) :: e_unit
  !> Fermi energy
  real(dp), intent(in) :: e_fermi
  !> Flattened and sorted scf-energies                                                                                                                                                                                                                                                                                                                             
  real(dp), intent(in) :: scf_energies(:)
  !> file name
  character(len=*), intent(in) :: file_name

  !local                                                                                                                                                                                                                                                                                                                                                          
  integer(i32) :: i_element
  integer(i32) :: file_unit

  open(newunit=file_unit, file=file_name)
  write(file_unit,'(A12)') "version: 100"
  write(file_unit,'(A12)') "type: Tensor"
  write(file_unit,'(A18)') "scalarType: Real64"
  write(file_unit,'(A11)') "dimensions:"
  write(file_unit,'(A9, 1X, I0)') "- length:", size(scf_energies)
  write(file_unit,'(2X, A11)') "type: State"
  write(file_unit,'(A9)') "elements:"
  write(file_unit,'(2X, A14)') "type: TextFile"
  write(file_unit,'(A5, 1X, F12.7)') "unit:", e_unit
  write(file_unit, '(A9)' ) "metaData:"
  write(file_unit,'(2X, A12, 1X, F25.16)') "fermiEnergy:", e_fermi
  write(file_unit,'(2X, A9)') "energies:"

  do i_element=1, size(scf_energies)
    write(file_unit,'(2X, A1, 1X, F25.16)') "-", scf_energies(i_element)
  end do

  close(unit=file_unit)

end subroutine write_eigenenergies_to_yaml


!> Writes eigenenergies to an output file. 
subroutine write_scf_energies_to_file(scf_energies, n_states, lowest_state, file_name)
  !> Number of states
  integer(i32), intent(in) :: n_states
  !> Lowest state
  integer(i32), intent(in) :: lowest_state
  !> Eigen energies
  real(dp), intent(in) :: scf_energies(n_states)
  !> File name
  character(len=*), intent(in) :: file_name

  !local variables
  integer(i32) :: i_state
  integer(i32) :: file_unit

  open(newunit=file_unit, file=trim(file_name), action='write')
    do i_state=lowest_state,n_states
      write(file_unit,'(F25.16)') scf_energies(i_state)
    end do
  close(unit=file_unit)

end subroutine write_scf_energies_to_file


!> Writes orbital properties to an output file. 
subroutine write_orbital_properties_to_yaml(ordered_indices, n_states, file_name)
  !> indices of size ordered eigenvalues
  integer(i32), dimension(:), intent(in) :: ordered_indices
  !> number of states
  integer(i32), intent(in) :: n_states
  !> file name
  character(len=*), intent(in) :: file_name


  !local
  integer(i32) :: i,j,lb
  integer(i32) :: file_unit
  integer(i32), allocatable :: unordered_k_indices(:)
  integer(i32), allocatable :: unordered_spin_indices(:)

  !parameters
  ! Make this input, when implemented for spin and k-points
  integer(i32), parameter :: n_k_points = 1, n_spins = 1 

  allocate(unordered_k_indices(n_k_points*n_states*n_spins))
  allocate(unordered_spin_indices(n_k_points*n_states*n_spins))

  !Momentum indices
  do i=1,n_k_points
    lb = (i-1)*n_states*n_spins + 1
    unordered_k_indices(lb:lb+n_states*n_spins-1) = i
  end do

  !Spin indices
  do i=1,n_k_points
    do j=1,n_spins
      lb = ((i-1)*n_spins + (j-1))*n_states + 1
      unordered_spin_indices(lb:lb+n_states-1) = j
    end do
  end do

  open(newunit=file_unit, file=file_name)
  write(file_unit,'(A12)') "version: 100"
  write(file_unit,'(A20)') "dimensionType: State"
  write(file_unit,'(A11)') "properties:"
  write(file_unit,'(2X, A5)') "Spin:"

  if(n_spins == 2) then
    write(file_unit,'(4X, A7)') "0: +0.5"
    write(file_unit,'(4X, A7)') "1: -0.5"
  else
    write(file_unit,'(4X, A7)') "0: +0.5"
  end if 

  write(file_unit,'(A16)') "propertyIndices:"
  write(file_unit,'(2X, A16)') "CrystalMomentum:"
  do i=1,size(ordered_indices)
    write(file_unit,'(2X, A1, 1X, I0)') "-", unordered_k_indices(ordered_indices(i)) - 1
  end do

  write(file_unit,'(2X, A5)') "Spin:"
  do i=1,size(ordered_indices)
    write(file_unit,'(2X, A1, 1X, I0)') "-", unordered_spin_indices(ordered_indices(i)) - 1
  end do

  close(file_unit)
end subroutine write_orbital_properties_to_yaml


!> Write a small test/reference file summarizing the Coulomb vertex for a
!> given q-point.
!>
!> The full vertex array is typically too large to store or diff directly
!> in a test suite, so this routine reduces it to cheap, deterministic
!> fingerprints instead: the array shape and its Frobenius norm
!> \f$ \sqrt{\sum_{ijk} |V_{ijk}|^2} \f$. Comparing these values against a
!> reference (with a numerical tolerance on the norm) is sufficient to catch
!> most regressions without needing to store the full array.
!>
!> The output is written to a plain-text file named
!> `coulomb_vertex_q<iq>.test`, with `iq` zero-padded to 4 digits.
subroutine write_test_output(iq, coulomb_vertex)
  !> current q-point index (right now, always 1)
  integer(i32), intent(in) :: iq
  !> Coulomb vertex
  complex(dp), intent(in)  :: coulomb_vertex(:, :, :)

  character(len=128) :: fname
  integer(i32) :: fid
  real(dp) :: checksum

  ! Frobenius norm of the vertex, used as a cheap fingerprint for testing
  checksum = sqrt(sum(abs(coulomb_vertex)**2))

  write(fname, '("coulomb_vertex_q", i4.4, ".test")') iq

  open(newunit=fid, file=trim(fname), status='replace', action='write')
  write(fid, '("iq        = ", i0)') iq
  write(fid, '("shape     = ", 3(i0, 1x))') shape(coulomb_vertex)
  write(fid, '("coulomb_vertex  = ", es23.16)') checksum
  close(fid)

end subroutine write_test_output


!> Compute the Coulomb vertex:
!> $\Gamma_r^{p\nu}=\sum_{\mu}M_r^{p\mu}\sqrt{v_\mu^\nu}$, 
!> with
!> $\psi_p^*(\mathbf{r})\psi_r(\mathbf{r})=\sum_{\mu}M_r^{p\mu}\chi_\mu(\mathbf{r})$
!> and
!> $v_\mu^\nu=\int\frac{\chi_\mu(\mathbf{r})\chi_\nu^*(\mathbf{r}')}
!> {|\mathbf{r}-\mathbf{r}'|}\,d\mathbf{r}\,d\mathbf{r}'$.
subroutine compute_coulomb_vertex(iq, iomstart, iomend, mbsiz, n_last_frozen)

  use constants, only: zone, zzero, pi
  use mod_APW_LO, only: apwordmax
  use mod_atoms, only: natmtot
  use mod_muffin_tin, only: lmmaxapw
  use modgw, only: kqset, fnm, fnm_tet, fnm_sum, mblksiz, b2mb, gkqset, msize
  use mod_product_basis, only: minmmat
  use mod_bands, only: nstdf, numin, eveckpalm, eveckalm, eveck, eveckp
  use mod_eigenvalue_occupancy, only: nstfv
  use mod_eigensystem, only: nmatmax 
  use mod_device_offload,    only: device_world
  use device_linalg_common_interface, only: zgemm_gpu
  use mod_selfenergy,        only: singc2
  use mod_expand_products, only: expand_products_generic
  use mod_misc_gw,           only: vi

  !> q-point index  
  integer(i32), intent(in) :: iq
  !> Initial and final frequency grid index
  integer(i32), intent(in) :: iomstart, iomend
  !> Mixed Product basis size
  integer(i32), intent(in) :: mbsiz
  !> Last index of the frozen states
  integer(i32), intent(in) :: n_last_frozen

  ! local
  integer(i32) :: ie1, ie2, ibasis
  integer(i32) :: ik, jk
  integer(i32) :: ndim
  integer(i32) :: nblk, iblk
  real(dp)    :: corr 
  complex(dp), allocatable :: evecfv(:,:)
  integer(i32) :: my_device
  logical :: tetrahedron_method

  !=============================
  ! Initialization
  !=============================

  ! Get device id
  my_device = device_world%get_device()

  ndim = nstdf - n_last_frozen + 1

  ! arrays to store products of KS eigenvectors with the matching coefficients
  allocate(eveckalm(nstfv,apwordmax,lmmaxapw,natmtot))
  allocate(eveckpalm(nstfv,apwordmax,lmmaxapw,natmtot))
  allocate(eveck(nmatmax,nstfv))
  allocate(eveckp(nmatmax,nstfv))

  OMP_OFFLOAD target enter data map(alloc: eveckalm, eveckpalm, eveck, eveckp)


  !=================
  ! BZ integration
  !=================

  ! k-q point
  ik = 1
  jk = kqset%kqid(ik, iq)


  ! get KS eigenvectors
  allocate(evecfv(nmatmax,nstfv))
  call get_evec_gw(kqset%vkl(:,jk), Gkqset%vgkl(:,:,:,jk), evecfv)
  eveckp = conjg(evecfv)
  call get_evec_gw(kqset%vkl(:,ik), Gkqset%vgkl(:,:,:,ik), evecfv)
  eveck = evecfv
  deallocate(evecfv)

  ! compute products \sum_G C_{k}n * A_{lm}
  call expand_evec(ik,'t')
  call expand_evec(jk,'c')

  OMP_OFFLOAD target update to(eveck, eveckp, eveckalm, eveckpalm)

  !=================================================
  ! Loop over m-blocks in M^i_{nm}(\vec{k},\vec{q})
  !=================================================

  allocate(minmmat(mbsiz, n_last_frozen:nstdf, n_last_frozen:nstdf))
  OMP_OFFLOAD target enter data map(alloc: minmmat)

  ! Compute M^i_{nm}+M^i_{cm}
  call expand_products_generic(ik, iq, n_last_frozen, nstdf, n_last_frozen, & 
                                0, n_last_frozen, nstdf, 1, 0, minmmat, .true.)

  ! ToDo: port to GPU
  coulomb_vertex(1:mbsiz, n_last_frozen:nstdf, n_last_frozen:nstdf) = minmmat

  deallocate(minmmat)

  ! To circumvent the divergence in v matrix elements, the constant G=0 component in 
  ! Chi_mu is discarded. To afterwards correct the Coulomb vertex, 
  ! the method by Massida et al. [Phys. Rev. B48, 5058 (1993)] is used. 
  ! This leads to following adjustment of the mu=1 component of the Coulomb vertex:
  ! $\Gamma_{p1}^{p}=\sqrt{\sqrt{\frac{1}{\pi\beta}}-\frac{4\pi}{\Omega}\sum_{G\neq 0}\frac{e^{-\beta G^2}}{G^2}}$ 
  corr = 4.0_dp * pi * vi * singc2
  do ie1 = 1, nstdf 
    coulomb_vertex(mbsiz+1, ie1, ie1) = cmplx(sqrt(corr), 0.0_dp)
  end do

  !-------------------
  ! Clear memory
  !-------------------

  OMP_OFFLOAD target update from(coulomb_vertex)

  OMP_OFFLOAD target exit data map(delete: eveck, eveckp, eveckalm, eveckpalm,minmmat)

  deallocate(eveck)
  deallocate(eveckp)
  deallocate(eveckalm)
  deallocate(eveckpalm)
end subroutine compute_coulomb_vertex


!> Prepare dft eigenvalues for cc4s.
subroutine prepare_dft_eigenvalues(nocc, nunocc, scf_eigenvalues, &
                                  scf_eigenvalues_flattened, &
                                  sorted_ids_scf_eigvals)

  use mod_bands, only: evalfv
  use sorting, only: sortidx

  !> Number of occupied / unoccpuied states
  integer(i32), intent(in) :: nocc, nunocc
  !> eigenvalues for every state, spin, and k-point
  real(dp), allocatable, intent(out) :: scf_eigenvalues(:, :, :)
  !> eigenvalues for every state, spin, and k-point
  !> flattend into one dimension
  real(dp), allocatable, intent(out) :: scf_eigenvalues_flattened(:)
  !> eigenvalues sorted by size
  integer(i32), allocatable, intent(out) :: sorted_ids_scf_eigvals(:)

  allocate(scf_eigenvalues_flattened(size(evalfv(:, 1))))
  ! Note: This only works for one spin and one k-point.
  !       Needs to be adjusted, when we can handle more.
  scf_eigenvalues_flattened(:) = evalfv(:, 1) 

  allocate(scf_eigenvalues(size(evalfv(:, 1)), 1, 1)) ! (n_states, n_spin, n_k)
  scf_eigenvalues(:, 1, 1) = evalfv(:, 1) ! (n_states, n_spin, n_k) = (n_states, n_spin)

  ! Get indices of sorted DFT eigenvalues
  allocate(sorted_ids_scf_eigvals(size(scf_eigenvalues_flattened)))
  call sortidx(size(scf_eigenvalues_flattened), scf_eigenvalues_flattened, sorted_ids_scf_eigvals)

end subroutine prepare_dft_eigenvalues


end module
