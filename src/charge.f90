!> Computes the muffin-tin, interstitial and total charges by integrating the
!> real-space density stored in arrays `rhomt` and `rhoir`.
subroutine charge( caller_tag )
  use constants, only: fourpi, y00
  use mod_atoms, only: nspecies, natoms, idxas, spr
  use mod_charge_and_moment, only: chgmttot, chgmt, chgtot, chgir, chgcalc
  use mod_Gvector, only: ngrtot, cfunir
  use mod_lattice, only: omega
  use mod_muffin_tin, only: nrmt, nrmtmax
  use mod_potential_and_density, only: rhomt, rhoir
  use modinput, only: input
  use precision, only: i32, dp
  use to_char_conversion, only: to_char

  implicit none

  !> String from the calling routine to be included in the warning message
  character(len=*), intent(in) :: caller_tag
  
  integer(i32) :: is, ia, ias
  real(dp) :: fr(nrmtmax), gr(nrmtmax), cf(3, nrmtmax) 

  chgmttot = 0._dp
  do is = 1, nspecies
    do ia = 1, natoms(is)
      ias = idxas(ia, is)
      fr(1 : nrmt(is)) = rhomt(1, 1 : nrmt(is), ias) * spr(1 : nrmt(is), is)**2
      call fderiv( -1, nrmt(is), spr(:, is), fr, gr, cf )
      chgmt(ias) = fourpi * y00 * gr(nrmt(is))
      chgmttot = chgmttot + chgmt(ias)
    end do
  end do
  
  chgir = dot_product( rhoir(1 : ngrtot), cfunir(1 : ngrtot) ) * omega / real( ngrtot, kind = dp )
  chgcalc = chgmttot + chgir

  if ( abs( chgtot / chgcalc - 1._dp ) > input%groundstate%epschg ) then
    call warning( 'Warning(charge): Charge ' // to_char( chgcalc ) // &
      ' calculated from electronic density in ' // trim( caller_tag ) // &
      ' differs from the total GS charge ' // to_char( chgtot ) )
  end if
end subroutine
