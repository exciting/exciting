!BOP
! 
! !ROUTINE: intern
!
! !INTERFACE:
      subroutine intern(nkp, kp, div, shift, divsh, klist, idiv)

! !DESCRIPTION:
!
! This subroutine transform the submesh coordinates of the kpoints into
! internal coordinates in the basis vectors of the reciprocal lattice
!
! !INPUT PARAMETERS:
!
      implicit none

      integer(4), intent(in) :: nkp       ! Number of k-points
      integer(4), intent(in) :: kp(3, nkp) ! submesh coordinates of the k-points
      integer(4), intent(in) :: div(3)    ! number of divisions of 
                                          ! the  submesh in each direction
      integer(4), intent(in) :: shift(3)  ! shift of the submesh from  the origin
      integer(4), intent(in) :: divsh     ! common divisor of the shift

! !OUTPUT PARAMETERS:

      integer(4), intent(out) :: klist(3, nkp)! integer coordinates of the kpoints
      integer(4), intent(out) :: idiv         ! minimum common divisor for
                                              ! the integer coordinates of the kpoints
      
! !LOCAL VARIABLES:

      integer(4) :: kpi
      integer(4) :: i
      real(8)    :: rind
      
!EOP
!
!BOC      

      ! Why 2?
      !idiv = div(1)*div(2)*div(3)*2
! Note: The i'th component of a k-vector is given by
!       k_i = (\tilde{k}_{i,sub} + shift_{i,sub})/(divsh*ngridk(i)),
!       where sub is refering to the finer submesh that includes the offset
!       To assing each k a integer index we scale by idiv:
      idiv = div(1)*div(2)*div(3)*divsh
      do kpi = 1, nkp
        do i = 1, 3
          ! Get the i'th coordinate of the reduced k point
          ! in lattice coordinates
! Note: These are the same as those produced by genppts
          !rind = dble(divsh*kp(i, kpi)+shift(i))/dble(divsh*div(i))
          !klist(i, kpi) = nint(rind*idiv)
! Note: Try scaled coordinates 
          rind=dble((divsh*kp(i, kpi)+shift(i))*idiv)/dble(div(i)*divsh)
          klist(i, kpi) = nint(rind)
        end do
      end do

      ! Reduce common divisior idiv to smallest common divisor
      call divisi(nkp, idiv, klist)

      end subroutine intern

!EOC
