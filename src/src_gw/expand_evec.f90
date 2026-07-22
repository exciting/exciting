!BOP
!
!!ROUTINE: \verb"expand_evec"
!
!!INTERFACE:
!
    subroutine expand_evec(ik,trans)
!
!!DESCRIPTION:
!
!Calculate the product of an eigenvector with the corresponding matching coefficients
!
!!USES:
    use mod_get_eigenvectors_times_matchingcoefficients, only: get_eigenvectors_times_matchingcoefficients
    use mod_bands, only: eveckalm, eveckpalm, eveck, eveckp

!!INPUT PARAMETERS:
    implicit none
    integer(4),   intent(in) :: ik    ! index of the k-point      
    character(1), intent(in) :: trans

    select case (trans)
    case ('t','T')
      call get_eigenvectors_times_matchingcoefficients(ik, trans, eveck, eveckalm)
    case ('c','C')
      call get_eigenvectors_times_matchingcoefficients(ik, trans, eveckp, eveckpalm)
    case default
      write(*,*)'ERROR in expand_evec'
      write(*,*)'Allowed values of trans are t, T, c or C'
      write(*,'(a19,a1)')'Received trans = ',trans 
      stop 
    end select
end subroutine
!EOC      
