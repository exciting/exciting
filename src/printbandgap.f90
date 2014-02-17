!
! Created Feb 2014 by DIN based on src_xs/writebandgap.f90
!
subroutine printbandgap(fid)
    use modmain
    use modmpi
    use m_bandgap
    implicit none
    ! input variables
    integer, intent(in) :: fid
    ! local variables
    integer :: ikgf(2), ikgo, istho
    real(8) :: egf, ego
    real(8) :: ebmax, ebmin, egap(3)
    logical :: metallic
    
    call bandgap(nkpt,evalsv,efermi,egf,ego,ikgf,ikgo,istho)
    
    write(fid,*)
    if (egf > 0.d0) then
      write(fid,'(a, T45,": ", F18.8)') ' Fundamental gap', egf
      write(fid,'(a, i5, 4x, 3f8.4, F18.8)') '  k-point (homo), energy: ', &
      &  ikgf(1), vkl(:,ikgf(1)), evalsv(istho,ikgf(1))
      write(fid,'(a, i5, 4x, 3f8.4, F18.8)') '  k-point (lumo), energy: ', &
      &  ikgf(2), vkl(:,istho+1), evalsv(istho+1,ikgf(2))
      write(fid,'(a, T45,": ", F18.8)') ' Optical gap', ego
      write(fid,'(a, i5, 4x, 3f8.4, F18.8)') '  k-point (homo), energy: ', &
      &  ikgo, vkl(:,ikgo), evalsv(istho,ikgo)
      write(fid,'(a, i5, 4x, 3f8.4, F18.8)') '  k-point (lumo), energy: ', &
      &  ikgo, vkl(:,ikgo), evalsv(istho+1,ikgo)
    end if

    return
end subroutine
