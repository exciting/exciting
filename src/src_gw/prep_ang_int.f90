!BOP
! !ROUTINE: prep_ang_int
!
! !INTERFACE:
      subroutine prep_ang_int(ml,ngrid)
!
! !DESCRIPTION:
!
! This subroutine sets the grids for the integration of spherical harmonic
! products using Lebedev-Laikov algorithm and calculates the spherical
! harmonics up to ml for those grid points.
! 
! !USES:

      use modgw, only: nleb, sphar, wleb

! !INPUT PARAMETERS:
!
      implicit none

      integer(4) :: ngrid ! estimated number of gridpoints
!
      integer(4) :: ml    ! maximum l for which the spherical harmonics
!                           are calculated
!
!
! !LOCAL VARIABLES:
      integer(4) :: dimsph ! dimension of the sphar vector (number of 
!                            spherical harmonics) = (ml+1^2)
!
      integer(4) :: ileb   ! index of the grid point
!
      real(8), allocatable :: xleb(:) ! temporary storage for the 
!                                       x coordinate of the gridpoints
!
      real(8), allocatable :: yleb(:) ! temporary storage for the 
!                                       y coordinate of the gridpoints
!
      real(8), allocatable :: zleb(:) ! temporary storage for the 
!                                       z coordinate of the gridpoints
!
      real(8), dimension(3) :: gp     ! stores one gridpoint
!      

      complex(8), allocatable :: sa(:) ! temporary storage for 
!                                                 the spherical harmonics
!
!EOP
!BOC
      dimsph = (ml+1)*(ml+1)
      allocate(sa(dimsph))

grid: select case (ngrid)
        case (:6)            
          allocate(xleb(6))
          allocate(yleb(6))
          allocate(zleb(6))
          call init_angint(dimsph,6)
          call ld0006(xleb,yleb,zleb,wleb,nleb)
        case (7:14)
          allocate(xleb(14))
          allocate(yleb(14))
          allocate(zleb(14))
          call init_angint(dimsph,14)
          call ld0014(xleb,yleb,zleb,wleb,nleb)
        case (15:26)
          allocate(xleb(26))
          allocate(yleb(26))
          allocate(zleb(26))
          call init_angint(dimsph,26)
          call ld0026(xleb,yleb,zleb,wleb,nleb)
        case (27:38)
          allocate(xleb(38))
          allocate(yleb(38))
          allocate(zleb(38))
          call init_angint(dimsph,38)
          call ld0038(xleb,yleb,zleb,wleb,nleb)
        case (39:50)
          allocate(xleb(50))
          allocate(yleb(50))
          allocate(zleb(50))
          call init_angint(dimsph,50)
          call ld0050(xleb,yleb,zleb,wleb,nleb)
        case (51:74)
          allocate(xleb(74))
          allocate(yleb(74))
          allocate(zleb(74))
          call init_angint(dimsph,74)
          call ld0074(xleb,yleb,zleb,wleb,nleb)
        case (75:86)
          allocate(xleb(86))
          allocate(yleb(86))
          allocate(zleb(86))
          call init_angint(dimsph,86)
          call ld0086(xleb,yleb,zleb,wleb,nleb)
        case (87:110)
          allocate(xleb(110))
          allocate(yleb(110))
          allocate(zleb(110))
          call init_angint(dimsph,110)
          call ld0110(xleb,yleb,zleb,wleb,nleb)
        case (111:146)
          allocate(xleb(146))
          allocate(yleb(146))
          allocate(zleb(146))
          call init_angint(dimsph,146)
          call ld0146(xleb,yleb,zleb,wleb,nleb)
        case (147:170)
          allocate(xleb(170))
          allocate(yleb(170))
          allocate(zleb(170))
          call init_angint(dimsph,170)
          call ld0170(xleb,yleb,zleb,wleb,nleb)
        case (171:194)
          allocate(xleb(194))
          allocate(yleb(194))
          allocate(zleb(194))
          call init_angint(dimsph,194)
          call ld0194(xleb,yleb,zleb,wleb,nleb)
        case (195:230)
          allocate(xleb(230))
          allocate(yleb(230))
          allocate(zleb(230))
          call init_angint(dimsph,230)
          call ld0230(xleb,yleb,zleb,wleb,nleb)
        case (231:266)
          allocate(xleb(266))
          allocate(yleb(266))
          allocate(zleb(266))
          call init_angint(dimsph,266)
          call ld0266(xleb,yleb,zleb,wleb,nleb)
        case (267:302)
          allocate(xleb(302))
          allocate(yleb(302))
          allocate(zleb(302))
          call init_angint(dimsph,302)
          call ld0302(xleb,yleb,zleb,wleb,nleb)
        case (303:350)
          allocate(xleb(350))
          allocate(yleb(350))
          allocate(zleb(350))
          call init_angint(dimsph,350)
          call ld0350(xleb,yleb,zleb,wleb,nleb)
        case (351:434)
          allocate(xleb(434))
          allocate(yleb(434))
          allocate(zleb(434))
          call init_angint(dimsph,434)
          call ld0434(xleb,yleb,zleb,wleb,nleb)
        case (435:590)
          allocate(xleb(590))
          allocate(yleb(590))
          allocate(zleb(590))
          call init_angint(dimsph,590)
          call ld0590(xleb,yleb,zleb,wleb,nleb)
        case (591:770)
          allocate(xleb(770))
          allocate(yleb(770))
          allocate(zleb(770))
          call init_angint(dimsph,770)
          call ld0770(xleb,yleb,zleb,wleb,nleb)
        case (771:974)
          allocate(xleb(974))
          allocate(yleb(974))
          allocate(zleb(974))
          call init_angint(dimsph,974)
          call ld0974(xleb,yleb,zleb,wleb,nleb)
        case (975:1202)
          allocate(xleb(1202))
          allocate(yleb(1202))
          allocate(zleb(1202))
          call init_angint(dimsph,1202)
          call ld1202(xleb,yleb,zleb,wleb,nleb)
        case (1203:1454)
          allocate(xleb(1454))
          allocate(yleb(1454))
          allocate(zleb(1454))
          call init_angint(dimsph,1454)
          call ld1454(xleb,yleb,zleb,wleb,nleb)
        case (1455:1730)
          allocate(xleb(1730))
          allocate(yleb(1730))
          allocate(zleb(1730))
          call init_angint(dimsph,1730)
          call ld1730(xleb,yleb,zleb,wleb,nleb)
        case (1731:2030)
          allocate(xleb(2030))
          allocate(yleb(2030))
          allocate(zleb(2030))
          call init_angint(dimsph,2030)
          call ld2030(xleb,yleb,zleb,wleb,nleb)
        case (2031:2354)
          allocate(xleb(2354))
          allocate(yleb(2354))
          allocate(zleb(2354))
          call init_angint(dimsph,2354)
          call ld2354(xleb,yleb,zleb,wleb,nleb)
        case (2355:2702)
          allocate(xleb(2702))
          allocate(yleb(2702))
          allocate(zleb(2702))
          call init_angint(dimsph,2702)
          call ld2702(xleb,yleb,zleb,wleb,nleb)
        case (2703:3074)
          allocate(xleb(3074))
          allocate(yleb(3074))
          allocate(zleb(3074))
          call init_angint(dimsph,3074)
          call ld3074(xleb,yleb,zleb,wleb,nleb)
        case (3075:3470)
          allocate(xleb(3470))
          allocate(yleb(3470))
          allocate(zleb(3470))
          call init_angint(dimsph,3470)
          call ld3470(xleb,yleb,zleb,wleb,nleb)
        case (3471:3890)
          allocate(xleb(3890))
          allocate(yleb(3890))
          allocate(zleb(3890))
          call init_angint(dimsph,3890)
          call ld3890(xleb,yleb,zleb,wleb,nleb)
        case (3891:4334)
          allocate(xleb(4334))
          allocate(yleb(4334))
          allocate(zleb(4334))
          call init_angint(dimsph,4334)
          call ld4334(xleb,yleb,zleb,wleb,nleb)
        case (4335:4802)
          allocate(xleb(4802))
          allocate(yleb(4802))
          allocate(zleb(4802))
          call init_angint(dimsph,4802)
          call ld4802(xleb,yleb,zleb,wleb,nleb)
        case (4803:5294)
          allocate(xleb(5294))
          allocate(yleb(5294))
          allocate(zleb(5294))
          call init_angint(dimsph,5294)
          call ld5294(xleb,yleb,zleb,wleb,nleb)
        case (5295:5810)
          allocate(xleb(5810))
          allocate(yleb(5810))
          allocate(zleb(5810))
          call init_angint(dimsph,5810)
          call ld5810(xleb,yleb,zleb,wleb,nleb)
        case (5811:)
          write(*,*)'rga: WARNING: Necessary number of points for exact &
     & Lebenev integration is greater than 5811'
          write(*,*)'ngrid = ',ngrid
          write(*,*)'The angular integrals are approximated'
          allocate(xleb(5810))
          allocate(yleb(5810))
          allocate(zleb(5810))
          call init_angint(dimsph,5810)
          call ld5810(xleb,yleb,zleb,wleb,nleb)
      end select grid
      do ileb=1,nleb
        gp(1)=xleb(ileb)
        gp(2)=yleb(ileb)
        gp(3)=zleb(ileb)
        call ylm(gp,ml,sa)
        sphar(ileb,1:dimsph)=sa(1:dimsph)
      enddo
      deallocate(xleb,yleb,zleb,sa)
      return

      contains
      
        subroutine init_angint(m,n)
        implicit none
          integer(4), intent(in) :: m,n
          allocate(wleb(n))
          allocate(sphar(n,m))
        end subroutine init_angint

        subroutine end_angint
          deallocate(wleb)
          deallocate(sphar)
        end subroutine end_angint

      end subroutine prep_ang_int
!EOC
