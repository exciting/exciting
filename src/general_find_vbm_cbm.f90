module general_find_vbm_cbm
    use precision, only: i32, dp
    use to_char_conversion, only: to_char
    use modmpi, only: terminate_if_false
    implicit none

    public :: find_vbm_cbm, &
              find_vbm_cbm_auxiliary

    interface find_vbm_cbm
        module procedure :: find_vbm_cbm_efermi, &
                            find_vbm_cbm_occupancy
    end interface find_vbm_cbm 

contains

    !> calculate vbm and cbm index with efermi
    subroutine find_vbm_cbm_efermi(start_band_ind,end_band_ind,nkpoint,e_band_KS,efermi,vbm_band_ind,cbm_band_ind,vbm_kpt_ind,cbm_kpt_ind,gap_min_kpt_ind)

        !> range of states
        integer(i32), intent(in)  :: start_band_ind, end_band_ind
        !> number of k-points          
        integer(i32), intent(in)  :: nkpoint      
        !> band energies        
        real(dp), intent(in)  :: e_band_KS(start_band_ind:,:) 
        !> Fermi energy
        real(dp), intent(in)  :: efermi   
        !> index of VBM
        integer(i32), intent(out) :: vbm_band_ind            
        !> index of CBM
        integer(i32), intent(out) :: cbm_band_ind            
        !> k-point index of VBM
        integer(i32), intent(out) :: vbm_kpt_ind            
        !> k-point index of CBM
        integer(i32), intent(out) :: cbm_kpt_ind            
        !> k-point index of min(VB-CB)
        integer(i32), intent(out) :: gap_min_kpt_ind            

        integer(i32) :: ik
        ! VBM band indices at each k point
        integer(i32) :: val_band_ind(nkpoint)
        ! CBM band indices at each k point
        integer(i32) :: cond_band_ind(nkpoint)
        ! VBM band energies at each k point
        real(dp) :: eho(nkpoint)
        ! CBM band energies at each k point
        real(dp) :: elu(nkpoint)
        
        !> search for VBM and CBM for each k-point
        do ik = 1, nkpoint
            val_band_ind(ik) = start_band_ind-1+count(e_band_KS(start_band_ind:end_band_ind,ik) <= efermi)
            cond_band_ind(ik) = val_band_ind(ik)+1
            eho(ik) = e_band_KS(val_band_ind(ik),ik)
            elu(ik) = e_band_KS(cond_band_ind(ik),ik)
        end do ! ik

        call find_vbm_cbm_auxiliary(start_band_ind, end_band_ind, val_band_ind, cond_band_ind, eho, elu, vbm_band_ind, cbm_band_ind, vbm_kpt_ind, cbm_kpt_ind, gap_min_kpt_ind)
        
    end subroutine

    !> calculate vbm and cbm index with occupancy
    subroutine find_vbm_cbm_occupancy(start_band_ind,end_band_ind,nkpoint,occ,e_band_KS,vbm_band_ind,cbm_band_ind,vbm_kpt_ind,cbm_kpt_ind,gap_min_kpt_ind)

        !> range of states
        integer(i32), intent(in)  :: start_band_ind, end_band_ind          
        !> number of k-points
        integer(i32), intent(in)  :: nkpoint              
        !> band occupancies 
        real(dp), intent(in)  :: occ(start_band_ind:,:)   
        !> band energies
        real(dp), intent(in)  :: e_band_KS(start_band_ind:,:) 
        !> index of VBM
        integer(i32), intent(out) :: vbm_band_ind            
        !> index of CBM
        integer(i32), intent(out) :: cbm_band_ind            
        !> k-point index of VBM
        integer(i32), intent(out) :: vbm_kpt_ind            
        !> k-point index of CBM
        integer(i32), intent(out) :: cbm_kpt_ind            
        !> k-point index of min(VB-CB)
        integer(i32), intent(out) :: gap_min_kpt_ind            

        integer(i32) :: i, ik

        ! VBM band indices at each k point
        integer(i32) :: val_band_ind(nkpoint)
        ! CBM band indices at each k point
        integer(i32) :: cond_band_ind(nkpoint)
        ! VBM band energies at each k point
        real(dp) :: eho(nkpoint)
        ! CBM band energies at each k point
        real(dp) :: elu(nkpoint)

        ! Threshold below which occupancy is treated as zero
        real(dp), parameter :: tolerance = 1.0e-8_dp
       
        ! search for VBM and CBM for each k-point
        do ik = 1, nkpoint
            val_band_ind(ik) = maxval(pack([(i, i=start_band_ind,end_band_ind)], occ(start_band_ind:end_band_ind,ik) >= tolerance))
            cond_band_ind(ik) = val_band_ind(ik)+1
            eho(ik) = e_band_KS(val_band_ind(ik),ik)
            elu(ik) = e_band_KS(cond_band_ind(ik),ik)
        end do ! ik

        call find_vbm_cbm_auxiliary(start_band_ind, end_band_ind, val_band_ind, cond_band_ind, eho, elu, vbm_band_ind, cbm_band_ind, vbm_kpt_ind, cbm_kpt_ind, gap_min_kpt_ind)

    end subroutine

    !> auxiliary subroutine to derive vbm and cbm indices
    subroutine find_vbm_cbm_auxiliary(start_band_ind, end_band_ind, val_band_ind, cond_band_ind, eho, elu, vbm_band_ind, cbm_band_ind, vbm_kpt_ind, cbm_kpt_ind, gap_min_kpt_ind)
        !> range of states
        integer(i32), intent(in)  :: start_band_ind, end_band_ind
        !> VBM band indices at each k point
        integer(i32), intent(in) :: val_band_ind(:)
        !> CBM band indices at each k point
        integer(i32), intent(in) :: cond_band_ind(:)
        !> VBM band energies at each k point
        real(dp), intent(in) :: eho(:)
        !> CBM band energies at each k point
        real(dp), intent(in) :: elu(:)

        !> index of VBM
        integer(i32), intent(out) :: vbm_band_ind
        !> index of CBM
        integer(i32), intent(out) :: cbm_band_ind
        !> k-point index of VBM
        integer(i32), intent(out) :: vbm_kpt_ind
        !> k-point index of CBM
        integer(i32), intent(out) :: cbm_kpt_ind
        !> k-point index of min(VB-CB)
        integer(i32), intent(out) :: gap_min_kpt_ind

        ! get VBM and CBM band indices
        vbm_band_ind = maxval(val_band_ind)
        cbm_band_ind = minval(cond_band_ind)

        ! get VBM and CBM k-point indices
        vbm_kpt_ind = maxloc(eho,dim=1)
        cbm_kpt_ind = minloc(elu,dim=1)
        gap_min_kpt_ind = minloc(elu-eho,dim=1)

        ! error control
        call terminate_if_false( vbm_band_ind >= start_band_ind, " ERROR(find_vbm_cbm): VBM is out of the specified band ranges! vbm_band_ind = " // to_char( vbm_band_ind ) // " < start_band_ind = " // to_char( start_band_ind ) )

        call terminate_if_false( vbm_band_ind <= end_band_ind, " ERROR(find_vbm_cbm): VBM is out of the specified band ranges! vbm_band_ind = " // to_char( vbm_band_ind ) // " > end_band_ind = " // to_char( end_band_ind ) )

        call terminate_if_false( cbm_band_ind >= start_band_ind, " ERROR(find_vbm_cbm): CBM is out of the specified band ranges! cbm_band_ind = " // to_char( cbm_band_ind ) // " < start_band_ind = " // to_char( start_band_ind ) )

        call terminate_if_false( cbm_band_ind <= end_band_ind, " ERROR(find_vbm_cbm): CBM is out of the specified band ranges! cbm_band_ind = " // to_char( cbm_band_ind ) // " > end_band_ind = " // to_char( end_band_ind ) )

    end subroutine

end module 
