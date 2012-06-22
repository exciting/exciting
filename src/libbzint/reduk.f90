
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: reduk
!
! !INTERFACE: 
subroutine reduk(nsymt,divsh,weight)
  !     
  ! !DESCRIPTION:
  !
  ! This subroutine calculates a set of reduced k-points for the integration
  ! with the tetrahedron method. The k-points are reduced by the symmetry and 
  ! a weight is put to each irreducible k-point to tell how many points it 
  ! represents before the symmetry operation.

  !
  ! !USES:

  use kgen_internals
  !<sag>
  use control, only: tetraifc, tetradbglv
  !</sag>

  implicit none
  !<sag>
  real(8), parameter :: epslat=1.d-6
  real(8) :: kpr(3),kprt(3)
  logical, allocatable :: done(:)
  !</sag>


  ! !INPUT PARAMETERS:

  integer(4), intent(in) :: nsymt     ! Number of symmetry operations
  integer(4), intent(in) :: divsh

  ! !OUTPUT PARAMETERS:

  integer(4), intent(out):: weight(*) ! weight of the 
  !                                           reducible k-point
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
  !      write(6,*)'------------------------------------------------------'
  !      write(6,*)'               reduk: begin'
  !      write(6,*)'------------------------------------------------------'

  !<sag>
  if (trim(tetraifc)=='wien2k') then

!!! *** OLD CODE WITH BUGS ***
!!$     ! original code
!!$     do i1=0,div(1)-1
!!$        kk(1)=divsh*i1+shift(1)
!!$        onek(1)=i1
!!$        do i2=0,div(2)-1
!!$           kk(2)=divsh*i2+shift(2)
!!$           onek(2)=i2
!!$           do i3=0,div(3)-1
!!$              kk(3)=divsh*i3+shift(3)
!!$              onek(3)=i3
!!$              kpid=idkp(onek)
!!$              if(jget(kpav,kpid).eq.0) then
!!$                 nirkp=nirkp+1
!!$                 ikpid(kpid)=nirkp
!!$                 do i=1,48
!!$                    starr(i)=0
!!$                 enddo
!!$                 !              write(24,*)
!!$                 !              write(24,'(i6,3i4)')divsh,div
!!$                 !              write(24,'(3i4,"  ",3i4)')onek,kk
!!$                 do i=1,nsymt
!!$                    do j=1,3
!!$                       nkp(j)=mod(iio(j,1,i)*kk(1)+iio(j,2,i)*kk(2)+&
!!$                            &                 iio(j,3,i)*kk(3),divsh*div(j))
!!$                       ktp(j)=nkp(j)
!!$                       nkp(j)=nkp(j)+(1-isign(1,nkp(j)))*divsh*div(j)/2
!!$                       nkp(j)=(nkp(j)-shift(j))/divsh
!!$                    enddo
!!$                    !               write(24,'(i6,3i4,"  ",3i4)')i,ktp,nkp
!!$                    nkpid=idkp(nkp)
!!$                    !                write(24,'(4x,2i4)')nirkp,nkpid
!!$                    call jset(kpav,nkpid,1)
!!$                    redkp(nkpid)=kpid
!!$                    starr(i)=nkpid
!!$                 enddo
!!$                 members=0
!!$                 do i=1,nsymt
!!$                    if(starr(i).ne.0) then
!!$                       members=members+1
!!$                       do j=i+1,nsymt
!!$                          if(starr(i).eq.starr(j)) starr(j)=0
!!$                       enddo
!!$                    endif
!!$                 enddo
!!$                 weight(nirkp)=members
!!$              else
!!$                 ikpid(kpid)=0
!!$              endif
!!$              ! find coords of reduced onek
!!$              call coorskp(redkp(kpid),nkp)
!!$           enddo
!!$        enddo
!!$     enddo

     allocate(done(div(1)*div(2)*div(3)))
     done(:)=.false.

     ! debug information
     if (tetradbglv.gt.0) then
        write(*,*) 'libbzint(reduk): number of symmetries: ',nsymt
        do i=1,nsymt
           write(*,*) i
           write(*,*) iio(1,:,i)
           write(*,*) iio(2,:,i)
           write(*,*) iio(3,:,i)
           write(*,*)
        end do
     end if

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
                 do i=1,48
                    starr(i)=0
                 enddo
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


  else if (trim(tetraifc)=='exciting') then
     ! new code

     allocate(done(div(1)*div(2)*div(3)))
     done(:)=.false.

     ! debug information
     if (tetradbglv.gt.0) then
        write(*,*) 'libbzint(reduk): number of symmetries: ',nsymt
        do i=1,nsymt
           write(*,*) i
           write(*,*) iio(1,:,i)
           write(*,*) iio(2,:,i)
           write(*,*) iio(3,:,i)
           write(*,*)
        end do
     end if

     do i3=0,div(3)-1
        kk(3)=divsh*i3+shift(3)
        onek(3)=i3
        do i2=0,div(2)-1
           kk(2)=divsh*i2+shift(2)
           onek(2)=i2
           do i1=0,div(1)-1
              kk(1)=divsh*i1+shift(1)
              onek(1)=i1
              kpid=idkp(onek)
              if (.not.done(kpid)) then
                 nirkp=nirkp+1
                 ikpid(kpid)=nirkp
                 do i=1,48
                    starr(i)=0
                 enddo
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

     ! end new code
  end if ! if (tetraifc)
  !</sag>

  ! *** DEBUG ***
  if (tetradbglv.gt.1) then
     write(*,*) 'REDUK REPORTS: kpid,ikpid(kpid),redkp(kpid)'
     write(*,*) 'div',div
     do kpid=1,div(1)*div(2)*div(3)
        write(*,*) kpid,ikpid(kpid),redkp(kpid)
     end do
     write(*,*)
  end if

  deallocate(iio)
  !      write(6,*)'------------------------------------------------------'
  !      write(6,*)'               reduk: end'
  !      write(6,*)'------------------------------------------------------'
  return

!<sag>
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
!</sag>
end subroutine reduk
!EOC
