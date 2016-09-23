!BOP
!
! !ROUTINE: reduk
!
! !INTERFACE: 
subroutine reduk(nsymt,divsh,weight)
     
! !DESCRIPTION:
! This subroutine calculates a set of reduced k-points for the integration
! with the tetrahedron method. The k-points are reduced by the symmetry and 
! a weight is put to each irreducible k-point to tell how many points it 
! represents before the symmetry operation.

! !USES:
  use kgen_internals

  implicit none
  real(8), parameter :: epslat=1.d-6
  real(8) :: kpr(3),kprt(3)
  logical, allocatable :: done(:)

  ! !INPUT PARAMETERS:
  integer(4), intent(in) :: nsymt     ! Number of symmetry operations
  integer(4), intent(in) :: divsh

  ! !OUTPUT PARAMETERS:
  integer(4), intent(out):: weight(*) ! weight of the reducible k-point

  ! !LOCAL VARIABLES:
  integer(4), allocatable :: kpav(:)

  integer(4) :: i,j,i1,i2,i3
  integer(4) :: onek(3),nkp(3),kk(3),kpid,nkpid,ktp(3)
  integer(4) :: members,starr(48)

  integer(4), external :: jget
  integer(4), external :: idkp

  !EOP
  !BOC
  external jset

  allocate(kpav(div(1)*div(2)*div(3)/32+1))

  nirkp=0
  do i=1,div(1)*div(2)*div(3)
     call jset(kpav,i,0)
  enddo

  allocate(done(div(1)*div(2)*div(3)))
  done(:)=.false.

  do i1=0,div(1)-1
     kk(1)=divsh*i1+shift(1)
     onek(1)=i1
     do i2=0,div(2)-1
        kk(2)=divsh*i2+shift(2)
        onek(2)=i2
        do i3=0,div(3)-1
           kk(3)=divsh*i3+shift(3)
           onek(3)=i3
           kpid=idkp(onek)
           if (.not.done(kpid)) then
              nirkp=nirkp+1
              ikpid(kpid)=nirkp
              starr(:)=0
              do i=1,nsymt
                 ! from internal coordinates to lattice coordinates
                 kpr(:)=dble(kk(:))/dble(divsh*div(:))
                 ! rotate k-point
                 kprt(:)=iio(:,1,i)*kpr(1)+iio(:,2,i)*kpr(2)+&
                      iio(:,3,i)*kpr(3)
                 ! map to reciprocal unit cell
                 call mapto01(epslat,kprt)
                 ! from lattice coordinates to internal coordinates
                 ! for unshifted mesh
                 kprt(:)=(kprt(:)*div(:)*divsh-shift(:))/dble(divsh)
                 ! location of rotated k-point on integer grid
                 nkp(:)=nint(kprt(:))
                 ! determine fractional part of rotated k-point
                 call mapto01(epslat,kprt)
                 ! if fractional part is present discard k-point
                 if (any(kprt.gt.epslat)) nkp(:)=-1
                 nkpid=idkp(nkp)
                 if (.not.all(nkp.eq.-1)) then
                    done(nkpid)=.true.
                    redkp(nkpid)=kpid
                    starr(i)=nkpid
                 end if
              enddo
              members=0
              do i=1,nsymt
                 if(starr(i).ne.0) then
                    members=members+1
                    do j=i+1,nsymt
                       if(starr(i).eq.starr(j)) starr(j)=0
                    enddo
                 endif
              enddo
              weight(nirkp)=members
           else
              ikpid(kpid)=0
           endif
           ! find coords of reduced onek
           call coorskp(redkp(kpid),nkp)
        enddo
     enddo
  enddo

  deallocate(done)
  deallocate(iio)

  return

contains

  subroutine mapto01(epsl,v)
    implicit none
    ! arguments
    real(8), intent(in) :: epsl
    real(8), intent(inout) :: v(3)
    ! local variables
    integer :: k
    do k=1,3
       v(k)=v(k)-dble(int(v(k)))
       if (v(k).lt.0.d0) v(k)=v(k)+1.d0
       if ((1.d0-v(k).lt.epsl).or.(v(k).lt.epsl)) v(k)=0.d0
    end do
  end subroutine mapto01

end subroutine reduk
!EOC
