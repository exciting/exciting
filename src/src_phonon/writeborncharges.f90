subroutine writeborncharges
  use modinput
  use mod_qpoint
  use mod_phonon
  use mod_atoms
  use mod_misc, only: filext

  integer :: iq, is, ia, ias, ip, j
  real(8) :: pol_tmp(3,3), pol(3,2,3,natmtot), zstar(3,3,0:natmtot), dph

  call init0
  ! reset file extension
  write( filext, '(".OUT")')
  ! find unperturbed polarization
  if( input%phonons%gamma == 'onestep') then
    call macro_polarization( pol_tmp, 'read')
    do ias = 1, natmtot
      do ip = 1, 3
        pol(:,1,ip,ias) = pol_tmp(:,1)
      end do
    end do
  end if
  ! find Gamma point and polarizations
  do iq = 1, nqpt
    if( all( ivq(:,iq) == 0)) then
      do is = 1, nspecies
        do ia = 1, natoms( is)
          ias = idxas( ia, is)
          do ip = 1, 3
            call phfext( iq, is, ia, ip, 0, 1, filext, filextdyn, phdirname)
            j = chdir( './'//trim( adjustl( phdirname)))
            call macro_polarization( pol_tmp, 'read')
            pol(:,2,ip,ias) = pol_tmp(:,1)
            if( input%phonons%gamma == 'onestep') then
              j = chdir( '..')
              cycle
            end if
            call phfext( iq, is, ia, ip, 0, 2, filext, filextdyn, phdirname)
            call macro_polarization( pol_tmp, 'read')
            pol(:,1,ip,ias) = pol_tmp(:,1)
            if( input%phonons%gamma == 'standard') pol(:,:,ip,ias) = -pol(:,:,ip,ias)
            j = chdir( '..')
          end do
        end do
      end do
      exit
    end if
  end do
  ! reset file extension
  write( filext, '(".OUT")')
  ! Gamma point not found
  if( iq > nqpt) return
  ! compute Born effective charges  
  dph = input%phonons%deltaph
  if( input%phonons%gamma .eq. 'twostep') dph = 2.d0*dph
  call borncharges( pol, dph, zstar, .true., 'calc')
  call borncharges( pol, dph, zstar, .true., 'write')

  return
end subroutine writeborncharges
