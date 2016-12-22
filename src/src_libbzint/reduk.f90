!BOP
!
! !ROUTINE: reduk
!
! !INTERFACE: 
subroutine reduk(nsymt, divsh, weight)
     
! !DESCRIPTION:
! This subroutine calculates a set of reduced k-points for the integration
! with the tetrahedron method. The k-points are reduced by the symmetry and 
! a weight is put to each irreducible k-point to tell how many points it 
! represents before the symmetry operation.

! !USES:
  use kgen_internals

  implicit none
  real(8), parameter :: epslat = 1.d-6
  real(8) :: kpr(3), kprt(3)
  logical, allocatable :: done(:)

  ! !INPUT PARAMETERS:
  integer(4), intent(in) :: nsymt     ! Number of symmetry operations
  integer(4), intent(in) :: divsh     ! Scaling factor of ngridk

  ! !OUTPUT PARAMETERS:
  integer(4), intent(out):: weight(*) ! weight of the reducible k-point

  ! !LOCAL VARIABLES:
  integer(4), allocatable :: kpav(:)

  integer(4) :: i, j, i1, i2, i3
  integer(4) :: onek(3), nkp(3), kk(3), kpid, nkpid, ktp(3)
  integer(4) :: members, starr(48)

  integer(4), external :: jget
  integer(4), external :: idkp

  !EOP
  !BOC
  external jset

! Note: kpav not used.
  ! One k-point repesents one bit in an 32bit integer array
  ! so one integer in the array contains 0/1 information about 
  ! a k point.
  allocate(kpav(div(1)*div(2)*div(3)/32+1))
  ! Set the first nkpt bits in kpav to 0
  do i = 1, div(1)*div(2)*div(3)
     call jset(kpav, i, 0)
  end do

  ! Zero counter for irreducible k points
  nirkp = 0

  ! Flag array 
  allocate(done(div(1)*div(2)*div(3)))
  done(:) = .false.

  ! Loop over 3d index of k points
  ! Note: 3rd dimension has fastest index, 1st slowest index
  ! Loop over k points on original grid
  do i1 = 0, div(1)-1
     ! Get corresponding point on the finer grid
     ! where the offset is integer
     kk(1) = divsh*i1+shift(1)
     ! Save index
     onek(1) = i1
     do i2 = 0, div(2)-1
        kk(2) = divsh*i2+shift(2)
        onek(2) = i2
        do i3 = 0, div(3)-1
           kk(3) = divsh*i3+shift(3)
           onek(3) = i3

           ! Get 1d index of original k-grid point.
           ! Note: The function relies on the above used order in the k loops
           kpid = idkp(onek)

           ! Check if the k-point is already represented in irreducible set
           if ( .not. done(kpid)) then

              nirkp = nirkp+1

              ! Set order number of irreducible k point
              ! i.e. Mapping between non-reduced k points and reduced ones
              !      where redundant k points are mapped to 0
              ikpid(kpid) = nirkp

              starr(:) = 0

              ! Loop over symmetry operations
              do i = 1, nsymt

                 ! from internal coordinates to lattice coordinates
                 ! i.e. from fine k-mensh to lattice coordinates
                 kpr(:) = dble(kk(:))/dble(divsh*div(:))

                 ! rotate k-point S*k=kr
                 kprt(:) = iio(:, 1, i)*kpr(1)&
                        &+ iio(:, 2, i)*kpr(2)&
                        &+ iio(:, 3, i)*kpr(3)

                 ! map to reciprocal unit cell [0,1)
                 call mapto01(epslat, kprt)

                 ! from lattice coordinates to internal coordinates
                 ! for unshifted mesh
! Note: This maps back to k-grid coordinates of the coarse grid 
                 kprt(:) = (kprt(:)*div(:)*divsh-shift(:))/dble(divsh)

                 ! location of rotated k-point on integer grid corresponding
                 ! to the fine k-mesh
                 nkp(:) = nint(kprt(:))

                 ! determine fractional part of rotated k-point
! Note: There should not be any, since the symmetries were
!       reduced to those which transform the offset onto a coarse grid
!       point.
                 call mapto01(epslat, kprt)

                 ! if fractional part is present discard k-point
                 if(any(kprt .gt. epslat)) nkp(:)=-1

                 ! Get 1d index of corresponding original k grid point
                 nkpid = idkp(nkp)

                 ! If everything went according to plan 
                 ! set the k-point corresponding to the rotated one 
                 ! to "done", write map ik->ikp and set starr value
                 if ( .not. all(nkp .eq. -1)) then
                    done(nkpid)=.true.
                    redkp(nkpid) = kpid
                    starr(i) = nkpid
                 end if

              ! Symmetry rotations on k point number kpid
              end do

              ! Check how many k-points where reached
              ! applying symmetry operations on k point number kpid
              members = 0
              do i = 1, nsymt
                 if(starr(i) .ne. 0) then
                    members = members+1
                    ! If one k point is reach using multiple sym ops:
                    ! Zero any duplicates, one transformation is all we need
                    do j = i+1, nsymt
                       if(starr(i) .eq. starr(j)) starr(j) = 0
                    end do
                 end if
              end do

              ! The weight of kpid
              weight(nirkp) = members

           ! k point is already represented in the previous
           else

              ! Map iknr -> ik = 0, i.e. not included in reduced set
              ikpid(kpid) = 0

           end if

! Note: Useless call ?
           ! find coords of reduced onek
           call coorskp(redkp(kpid), nkp)

        ! k grid loops
        end do
     end do
  end do

  deallocate(done)
  deallocate(iio)

  return

contains

  subroutine mapto01(epsl, v)
    implicit none
    ! arguments
    real(8), intent(in) :: epsl
    real(8), intent(inout) :: v(3)
    ! local variables
    integer :: k
    do k = 1, 3
       v(k) = v(k)-dble(int(v(k)))
       if (v(k).lt.0.d0) v(k) = v(k)+1.d0
       if ((1.d0-v(k).lt.epsl).or.(v(k).lt.epsl)) v(k) = 0.d0
    end do
  end subroutine mapto01

end subroutine reduk
!EOC
