
subroutine find_vbm_cbm(ib,nb,nk,eband,efermi,ibvm,ibcm,ikvm,ikcm,ikvc)

    implicit none
    integer, intent(in)  :: ib, nb          ! range of states
    integer, intent(in)  :: nk              ! number of k-points
    real(8), intent(in)  :: eband(ib:nb,nk) ! band energies
    real(8), intent(in)  :: efermi          ! Fermi energy
    integer, intent(out) :: ibvm            ! index of VBM
    integer, intent(out) :: ibcm            ! index of CBM
    integer, intent(out) :: ikvm            ! k-point index of VBM
    integer, intent(out) :: ikcm            ! k-point index of CBM
    integer, intent(out) :: ikvc            ! k-point index of min(VB-CB)
   
    ! local variables
    integer :: i, ik
    integer :: ibv(nk), ibc(nk)
    real(8) :: eho(nk), elu(nk)
    integer :: ikvm0(1), ikcm0(1), ikvc0(1)
    
    ! search for VBM and CBM for each k-point
    do ik = 1, nk
      ibv(ik) = ib-1+count(eband(ib:nb,ik) <= efermi)
      ibc(ik) = ibv(ik)+1
      eho(ik) = eband(ibv(ik),ik)
      elu(ik) = eband(ibc(ik),ik)
    end do ! ik
    
    ! VB maximum
    ibvm = maxval(ibv)
    ! CB minimum
    ibcm = minval(ibc)
    
    ikvm0 = maxloc(eho)
    ikcm0 = minloc(elu)
    ikvc0 = minloc(elu-eho)
    
    ikvm = ikvm0(1)
    ikcm = ikcm0(1)
    ikvc = ikvc0(1)
    
    ! error control
    if (ibvm < ib) then
        write(*,*) "ERROR(find_vbm_cbm): VBM is out of the specified band ranges!"
        write(*,*) "ibvm = ", ibvm, " < ib", ib
        stop
    endif
    if (ibvm > nb) then
        write(*,*) "ERROR(find_vbm_cbm): VBM is out of the specified band ranges!"
        write(*,*) "ibvm = ", ibvm, " > nb = ", nb
        stop
    endif
    if (ibcm < ib) then
        write(*,*) "ERROR(find_vbm_cbm): CBM is out of the specified band ranges!"
        write(*,*) "ibcm = ", ibcm, " < ib", ib
        stop
    endif
    if (ibcm > nb) then
        write(*,*) "ERROR(find_vbm_cbm): CBM is out of the specified band ranges!"
        write(*,*) "ibcm = ", ibcm, " > nb = ", nb
        stop
    endif
    
    return
end subroutine
