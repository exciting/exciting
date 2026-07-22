
!--------------------------------!
! DFT groundstate related stuff  !
!--------------------------------!

module mod_bands
#include "asserts.fpp"
    use gw_info, only: write_to_gwinfo, write_to_gwinfo_boxmessage
    use modgw, only : kset, fgw, hatree_2_eV => hev
    use modmpi, only : rank, terminate
    use precision, only: i32, dp, str_64
    use to_char_conversion, only: to_char
    use general_find_vbm_cbm, only: find_vbm_cbm

    implicit none

    private
    public :: evalfv, &
              occfv, &
              eveck, &
              eveckp, &
              eveckalm, &
              eveckpalm, &
              nomax, &
              ikvbm, &
              numin, &
              ikcbm, &
              ikvcm, &
              nstdf, &
              nstse, &
              metallic, &
              nkp1, &
              kvecs1, &
              eks1, &
              eqp1, &
              nkp2, &
              kvecs2, &
              eks2, &
              eqp2, &
              delete_bands, &
              bandstructure_analysis

! First-variational eigenvalues
    real(dp), allocatable :: evalfv(:,:)

! First-variational occupations
    real(dp), allocatable :: occfv(:,:)

! Eigenvectors at k
    complex(dp), allocatable :: eveck(:,:)
      
! Eigenvectors at k'=k-q
    complex(dp), allocatable :: eveckp(:,:)
      
! Spherical harmonic expansion coefficients at k
    complex(dp), allocatable :: eveckalm(:,:,:,:)
      
! Spherical harmonic expansion coefficients at k'=k-q
    complex(dp), allocatable :: eveckpalm(:,:,:,:)
    
! Position of Valence Band Maximum (VBM)      
    integer(i32) :: nomax
    integer(i32) :: ikvbm
    
! Position of Conduction Band Minimum (CBM)      
    integer(i32) :: numin
    integer(i32) :: ikcbm
    
! Position of the direct v->c (optical) gap
    integer(i32) :: ikvcm

! Number of states used to calculate the dielectric function
    integer(i32) :: nstdf

! Number of states used to calculate the self-energy
    integer(i32) :: nstse
    
! Metallicity flag
    logical :: metallic

!---------------------------------------------------------------
! To be used in the interpolation routine (band structure plot) 
!---------------------------------------------------------------

! Input 
    integer(i32) :: nkp1
    real(dp), allocatable :: kvecs1(:,:)
    real(dp), allocatable :: eks1(:,:), eqp1(:,:)

! Output (interpolated)
    integer(i32) :: nkp2
    real(dp), allocatable :: kvecs2(:,:)
    real(dp), allocatable :: eks2(:,:), eqp2(:,:)
    
contains

    subroutine delete_bands
        if (allocated(eveck)) deallocate(eveck)
        if (allocated(eveckp)) deallocate(eveckp)
        if (allocated(eveckalm)) deallocate(eveckalm)
        if (allocated(eveckpalm)) deallocate(eveckpalm)
    end subroutine

    !> Obtain the band gap, check if it is direct or indirect, and print the results to GW_INFO.OUT
    subroutine bandstructure_analysis(title, first_band, eigs, e_fermi, obtain_dos, kpt_indexes, kpt_lattice_coord)
      !> String to be printed to GW_INFO
      character(len=*), intent(in) :: title
      !> Index of the first band
      integer(i32), intent(in) :: first_band
      !> Set of KS or QP eigenvalues
      real(dp), intent(in) :: eigs(first_band:, :)
      !> Estimate of the Fermi energy
      real(dp), intent(in) :: e_fermi
      !> If `.true.`, then obtain DOS at the fermi level
      logical, intent(in) :: obtain_dos
      !> Indexes of k-points
      integer(i32), intent(in), optional :: kpt_indexes(:)
      !> Lattice coordinates of the k-points. It must be present if `kpt_indexes` is present
      real(dp), intent(in), optional :: kpt_lattice_coord(:, :)
      
      integer(i32) :: n_kpt, last_band, i_VBM, i_CBm
      real(dp) :: largest_eig, smallest_eig, e_gap
      real(dp) :: dos_fermi, kpt_VBM(3), kpt_CBm(3)
  
      real(dp), external :: dostet_exciting

      n_kpt = size( eigs, 2 )
      last_band = ubound( eigs, 1 )

      if( present(kpt_indexes) ) then 
        CALL_ASSERT( size(kpt_indexes) == n_kpt, 'kpt_indexes must have nkpt elements' )
        CALL_ASSERT( present(kpt_lattice_coord), 'kpt_lattice_coord must be present if kpt_indexes is present' ) 
        CALL_ASSERT( size(kpt_lattice_coord, 1)==3, 'kpt_lattice_coord must have size 3 along 1st dim.' )
        CALL_ASSERT( size(kpt_lattice_coord, 2)==n_kpt, 'kpt_lattice_coord must have size nkpt along 2nd dim.' )
        CALL_ASSERT( .not. obtain_dos, 'test_dos must be .false. if kpt_indexes is present')
      end if
  
      ! check Fermi energy for correspondence to the specified band range
      smallest_eig = minval(eigs)
      largest_eig = maxval(eigs)
      if ((largest_eig < e_fermi) .or. (smallest_eig > e_fermi)) then
        call terminate( "ERROR(bandstructure_analysis): Fermi energy is outside the specified electronic bands energy range!" )
      end if
  
      ! Search for the indices of VBM and CBM
      call find_vbm_cbm(first_band, last_band, n_kpt, eigs, e_fermi, nomax, numin, i_VBM, i_CBm, ikvcm)
      if( present(kpt_indexes) ) then
        ikvbm = kpt_indexes(i_VBM); ikcbm = kpt_indexes(i_CBm)
        kpt_VBM(:) = kpt_lattice_coord(:, i_VBM); kpt_CBm = kpt_lattice_coord(:, i_CBm)
      else
        ikvbm = i_VBM; ikcbm = i_CBm
        kpt_VBM(:) = kset%vkl(:, i_VBM); kpt_CBm = kset%vkl(:, i_CBm)
      end if
  
      ! Calculate DOS at the fermi level
      if( obtain_dos ) then
        dos_fermi = dostet_exciting(last_band-first_band+1, n_kpt, eigs, kset%ntet,kset%tnodes,kset%wtet,kset%tvol, e_fermi)
      else
        dos_fermi = 0._dp
      end if
  
      ! check for VBM and CBM overlap (metal)
      metallic = ( (nomax >= numin) .or. (abs(dos_fermi)>1.d-4) )
  
      if (rank==0) then
        call write_to_gwinfo_boxmessage( '-', trim(title) )
        call wrapper_write_gw_info( "Fermi energy: ", "(A,F10.4)", [e_fermi] )
        call wrapper_write_gw_info( "Energy range: ", "(A,2F10.4)", [smallest_eig, largest_eig] )
        call write_to_gwinfo( "Band index of VBM: " // to_char( nomax ) )
        call write_to_gwinfo( "Band index of CBm: " // to_char( numin ) )
        call write_to_gwinfo('')
        if (metallic) then
          call wrapper_write_gw_info("DOS at Fermi level: ", "(A,F8.4)", [dos_fermi])
          call write_to_gwinfo( "WARNING(bandstructure_analysis): Valence and Conduction bands overlap (metal)!" )
        else
          e_gap = eigs(numin, i_CBm) - eigs(nomax, i_VBM)
          if (i_CBm == i_VBM) then
            call wrapper_write_gw_info( "Direct BandGap (eV):", "(A,T40,F10.4)", [e_gap*hatree_2_eV])
            call wrapper_write_gw_info( "at k      = ", "(A, 3F8.3,' ik = " // to_char(ikvbm) // "')", kpt_VBM )
          else
            call wrapper_write_gw_info( "Indirect BandGap (eV):", "(A,T40,F10.4)", [e_gap*hatree_2_eV])
            call wrapper_write_gw_info( "at k(VBM) = ", "(A, 3F8.3,' ik = " // to_char(ikvbm) // "')", kpt_VBM )
            call wrapper_write_gw_info( "   k(CBm) = ", "(A, 3F8.3,' ik = " // to_char(ikcbm) // "')", kpt_CBm )
            call wrapper_write_gw_info( "Direct Bandgap at k(VBM) (eV):", "(A,T40,F10.4)", [(eigs(numin, i_VBM) - eigs(nomax, i_VBM))*hatree_2_eV] )
            call wrapper_write_gw_info( "Direct Bandgap at k(CBm) (eV):", "(A,T40,F10.4)", [(eigs(numin, i_CBm) - eigs(nomax, i_CBm))*hatree_2_eV] )
          end if
        end if
        call flushifc(fgw)
      end if
    end subroutine

    !> (private)
    subroutine wrapper_write_gw_info( message_beginning, string_format, x )
      character(len=*), intent(in) :: message_beginning
      character(len=*), intent(in) :: string_format
      real(dp), intent(in) :: x(:)

      character(len=str_64) :: string

      write( string, string_format ) message_beginning, x
      call write_to_gwinfo( trim(string) )
    end subroutine
  
    
end module
