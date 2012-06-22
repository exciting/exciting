
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: sym2int
!
! !INTERFACE:
      subroutine sym2int(rbas,nsym,symopc,symopi)
!
! !DESCRIPTION:
!
!    Redefines symmetry operation from cartesian coordinates to internal
!    coordinates of the Bravais lattice with a
!    unitary transformation  $\bar{S}' =\bar{u}^{-1} \bar{s} \bar{u}$
!
! !INPUT PARAMETERS:

      implicit none
      
      integer(4), intent(in) :: nsym          ! Number of symmetry  operations
      integer(4), intent(in) :: symopc(3,3,*) ! The symetry operations in cartesian coordinates
      real(8),    intent(in) :: rbas(3,3)     ! The basis vectors of the Bravais lattice

! !OUTPUT PARAMETERS:

      integer(4), intent(out) :: symopi(3,3,*) ! The symetry operations in internal coordinates

! !LOCAL VARIABLES:

      integer(4) :: i
      integer(4) :: ind
      integer(4) :: j

      real(8), dimension(3,3) :: a
      real(8), dimension(3,3) :: b
      real(8), dimension(3,3) :: gbas
      real(8), dimension(3,3) :: gbas1
      
! !DEFINED PARAMETERS:

      real(8), parameter :: pi = 3.14159265358979323846      

! !SYSTEM ROUTINES:
      
      intrinsic matmul
      
      intrinsic nint
      
      external gbass
!
!EOP
!
!BOC
!
!     Check that the number of symmetries does not exceed 48
!
      if(nsym.gt.48) then
       write(*,*) 'nsym  .gt. 48)',nsym
       stop 'sdefl nsym > 48'
      endif
!
!     Calculate the reciprocal lattice vectors
!
      call gbass(rbas,gbas)
!
!     Transpose gbas and divide by 2pi.
!
      do i=1,3
        do j=1,3
          gbas1(j,i)=gbas(i,j)/2.d0/pi   ! since gbass does not give the factor of 1/(2*Pi)
        enddo
      enddo

      do ind=1,nsym
        do i=1,3
          do j=1,3
            a(i,j)=symopc(i,j,ind)
          enddo
        enddo
        b=matmul(rbas,a)
        a=matmul(b,gbas1)
        do i=1,3
          do j=1,3
            symopi(i,j,ind)=nint(a(i,j))
          enddo
        enddo
      enddo

      end subroutine sym2int
!EOC
