!> This module declares the R-vectors and stars for the Fourier transform
module mod_fourier_interpolation_domain

    use precision, only: i32, dp

    private
    public :: init_real_space_domain_for_fourier_interpolation

    !> Number of R vectors
    integer(i32), public, protected :: nrr
    !> Integer indexes for R
    integer(i32), public, protected, allocatable :: rindex(:,:)
    !> Integer indexes for R considering symmetry
    integer(i32), public, protected, allocatable :: Sigma_rindex(:,:,:)

    !> real-space cutoff parameter (maximal length)
    real(dp) :: rmax

    !> basis vectors
    real(dp), public, protected :: rbas(3,3)

    !> Tolerance for zero
    real(dp), parameter :: zero_tolerance = 1.0e-6_dp

    ! Flag to check if the module variables have already initialized
    logical(i32) :: initialized_domain = .false.
      
contains

  !> Sets the indexes of the real space lattice vectors used for
  !> Fourier interpolation ordered by increasing length
  subroutine init_real_space_domain_for_fourier_interpolation

      use modinput
      use modmain

      implicit none

      integer(i32) :: i         ! (Counter): runs over coordinates
      integer(i32) :: ippw      ! (Counter): runs over plane waves
      integer(i32) :: ir1       ! (Counter): run over x-coord of G
      integer(i32) :: ir2       ! (Counter): run over y-coord of G
      integer(i32) :: ir3       ! (Counter): run over z-coord of G
      integer(i32) :: jppw      ! (Counter): runs over plane waves
      integer(i32) :: nr1       ! Max. ir1
      integer(i32) :: nr2       ! Max. ir2
      integer(i32) :: nr3       ! Max. ir3
      integer(i32) :: nr        ! Maximum numb,er of plane waves in intipw
      integer(i32) :: isym, npoint_group, ipoint_group
      integer(i32) :: irvec(3)  ! Integer coordinates of the G-vector
      integer(i32), allocatable  :: invrindex(:,:,:)

      real(dp)                  :: rr
      real(dp)                  :: rvec(3)   ! Cartesian coordinates of the R-vector
      integer(i32), allocatable :: rind(:,:) ! Temporary storage for rindex
      real(dp), allocatable     :: rlen(:)   ! Temporary storage for all the  qpg's

      integer(i32) :: ierr
      character(len=*), parameter :: sname="setrindex"

      ! If already done do nothing
      if (initialized_domain) return

      ! shortcut for basis vectors
      rbas(:,1) = input%structure%crystal%basevect(:,1)
      rbas(:,2) = input%structure%crystal%basevect(:,2)
      rbas(:,3) = input%structure%crystal%basevect(:,3)

      ! real-space cutoff parameter (maximal length)
      if (.not.associated(input%gw)) then
        rmax = 120.d0
      else
        rmax = input%gw%rmax
      end if

      rvec(:) = norm2(rbas(:,:),dim=1)

      nr1 = 2*nint(rmax/rvec(1))
      nr2 = 2*nint(rmax/rvec(2))
      nr3 = 2*nint(rmax/rvec(3))
      nr = (2*nr1+1)*(2*nr2+1)*(2*nr3+1)

      !--------- R vectors generation ----------

      allocate(rlen(1:nr),rind(3,1:nr))

      ippw = 0
      do ir1 = -nr1, nr1
        irvec(1) = ir1
        do ir2 = -nr2, nr2
          irvec(2) = ir2
          do ir3 = -nr3, nr3
            irvec(3) = ir3
            ! Transform to cartesian coordinates
            rvec = matmul(rbas,irvec)
            rr = norm2(rvec)
            if (rr<=rmax) then
              ippw = ippw + 1
              rlen(ippw) = rr
              rind(1:3,ippw) = irvec(1:3)
            end if
          end do
        end do
      end do
      nrr = ippw

      ! sort by increasing length using shell algorithm
      call shelsort(nrr,rind,rlen)

      allocate(rindex(3,1:nrr), source=rind(:,1:nrr))

      initialized_domain = .true.

      npoint_group = count(norm2(vtlsymc(:,1:nsymcrys),dim=1) < zero_tolerance)
      allocate(Sigma_rindex(3, npoint_group, 1:nrr), source=0)

      do ippw = 1, nrr
        ipoint_group = 0
        do isym = 1, nsymcrys
          if (norm2(vtlsymc(:,1:isym)) >= zero_tolerance) cycle
          ipoint_group = ipoint_group + 1
          Sigma_rindex(:, ipoint_group, ippw) = matmul(symlat(:,:,lsplsymc(isym)),rindex(:,ippw))
        end do
      end do

      deallocate(invrindex,rlen,rind)

  end subroutine init_real_space_domain_for_fourier_interpolation

end module mod_fourier_interpolation_domain

