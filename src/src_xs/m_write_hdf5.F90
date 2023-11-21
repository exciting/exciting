! Copyright(C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_write_hdf5
  use xhdf5, only: xhdf5_type
  use os_utils, only: join_paths

  implicit none

  contains

    subroutine write_spectra_hdf5(iq, foff, w, eps, loss, sigma, gname)
      use mod_lattice, only: omega
      use mod_misc, only: version
      use modinput
      use modmpi
      use modbse, only: nk_bse
      use modxs, only: escale, unitout, unit1
      use modxs, only: ivgmt, vqlmt, vgcmt, vqcmt
      use modxs, only: sptclg, ivgigq
      use fox_wxml
      use m_getunit
      use mod_hdf5
      use os_utils

      implicit none

      ! Arguments
      integer, intent(in) :: iq
      logical, intent(in) :: foff
      real(8), intent(in) :: w(:)
      complex(8), intent(in) :: eps(:,:,:)
      real(8), intent(in) :: loss(:,:,:)
      complex(8), intent(in) :: sigma(:)
      character(128), intent(in) :: gname 

      ! Local variables
      type(xmlf_t), save :: xf
      character(256) :: buffer
      character(*), parameter :: thisnam = 'writeeps'
      integer :: i, iw, igqmt
      real(8) :: w_(size(w))
      character(:), allocatable :: gname_, group, epsname
      character(4) :: momentum_index, ci
      complex(8), allocatable :: eps_(:)
      real(8), allocatable :: loss_(:)

      type(xhdf5_type) :: h5 


      !Call kramkron(iop1, iop2, 1.d-8, n, w, imeps, kkeps)

      igqmt = ivgigq(ivgmt(1,iq),ivgmt(2,iq),ivgmt(3,iq),iq)
      gname_ = gname

      call h5%initialize(fhdf5, mpiglobal, serial_access=.true.)
      call h5%initialize_group('.', gname_)

        write(momentum_index, '(I4.4)') iq ! Generate string out of momentum transfer index
        call h5%initialize_group(gname_, momentum_index)
        gname_ = join_paths(gname_, momentum_index)

        call h5%initialize_group(gname_, 'parameters')
        group = join_paths(gname_, 'parameters')

        call h5%write(group, "ivgmt", ivgmt(:, iq), [1], [3])
        call h5%write(group, "vqlmt", vqlmt(:, iq), [1], [3])
        call h5%write(group, "vgcmt", vgcmt(:, iq), [1], [3])
        call h5%write(group, "vqcmt", vqcmt(:, iq), [1], [3])
        call h5%write(group, "escale", escale)
        call h5%write(group, "broad", escale * input%xs%broad)
        call h5%write(group, "nk_bse", nk_bse)
        call h5%write(group, "foff", foff)

        call h5%initialize_group(gname_, 'diel')
        group = join_paths(gname_, 'diel')

        w_ = w * escale

        call h5%write(group, "w", w_, [1], shape(w_))

        if (foff) then ! Write full epsilon
          call h5%write(group, "epsm", eps, [1, 1, 1], shape(eps))
        else
          do i=1, 3 ! Write diagonal of epsilon
            write(ci, '(I4.2)') i*10+i
            epsname='epsm('//trim(adjustl(ci))//')'
            eps_ = eps(i, i, :)
            call h5%write(group, epsname, eps_, [1], shape(eps_))
          end do 
        end if 

        call h5%initialize_group(gname_, 'loss')
        group = join_paths(gname_, 'loss')

        call h5%write(group, "w", w_, [1], shape(w_))

        if (foff) then ! Write full loss function
          call h5%write(group, "lossfct", loss, [1, 1, 1], shape(loss))
        else
          do i=1, 3 ! Write diagonal of epsilon
            write(ci, '(I4.2)') i*10+i
            epsname='lossfct('//trim(adjustl(ci))//')' 
            loss_ = loss(i, i, :)
            call h5%write(group, epsname, loss_, [1], shape(loss_))
          end do 
        end if 
  
        call h5%initialize_group(gname_, 'sigma')
        group = join_paths(gname_, 'sigma')
        
        call h5%write(group, "w", w_, [1], shape(w_))
        call h5%write(group, "sigma", sigma, [1], shape(sigma))

        call h5%finalize()
   end subroutine write_spectra_hdf5

   subroutine write_excitons_hdf5(hamsize, nexc, eshift, evalre, oscstrr,&
      & gname, iqmt)
      use mod_lattice, only: omega
      use mod_misc, only: version
      use modinput
      use modmpi
      use modbse, only: nk_bse
      use modxs, only: escale, unitout, unit1
      use modxs, only: ivgmt, vqlmt, vgcmt, vqcmt
      use modxs, only: sptclg, ivgigq
      use fox_wxml
      use m_getunit
      use mod_hdf5
      use xhdf5, only: xhdf5_type
      use os_utils, only: join_paths

      implicit none

      ! Arguments
      integer(4), intent(in) :: hamsize, nexc
      real(8), intent(in) :: eshift
      real(8), intent(in) :: evalre(hamsize)
      complex(8), intent(in) :: oscstrr(:, :, :)
      character(128), intent(in) :: gname
      integer(4), intent(in) :: iqmt
   
      ! Local
      logical :: fsort
      integer(4) :: o1, o2, lambda, unexc, i, io1, io2, io3, io4, iq
      integer(4), allocatable :: idxsort(:), idxsort2(:)
      real(8), allocatable :: evalre_sorted(:)
      real(8), allocatable :: evalre_(:)
      
      character(:), allocatable :: gname_, group
      character(4) :: momentum_index, ci
      
      type(xhdf5_type) :: h5 

      gname_ = gname

      call h5%initialize(fhdf5, mpiglobal, serial_access=.true.)
      call h5%initialize_group('.', gname_)
      
      write(momentum_index, '(I4.4)') iqmt
      call h5%initialize_group(gname_, momentum_index)
      gname_ = join_paths(gname_, momentum_index)

      call h5%initialize_group(gname_, 'parameters')
      group = join_paths(gname_, 'parameters')

      iq = iqmt 

      call h5%write(group, "ivgmt", ivgmt(:, iq), [1], [3])
      call h5%write(group, "vqlmt", vqlmt(:, iq), [1], [3])
      call h5%write(group, "vgcmt", vgcmt(:, iq), [1], [3])
      call h5%write(group, "vqcmt", vqcmt(:, iq), [1], [3])
      call h5%write(group, "escale", escale)
      call h5%write(group, "eshift", eshift*escale)

      evalre_ = evalre * escale
      call h5%write(gname_, "evalre", evalre_, [1], shape(evalre_))
      deallocate(evalre_)

      io1=1
      io2=3
      if(iq /= 1) then 
        io1=1
        io2=1
      end if
      
      do o1 = io1, io2
        if ((iq .eq. 1) .and. (input%xs%dfoffdiag)) then
          io3=1
          io4=3
        else
          io3=o1
          io4=o1
        end if
        do o2 = io3, io4
          write(ci, '(I1.1, I1.1)') o1, o2
          call h5%initialize_group(gname_, trim(adjustl(ci)))
          call h5%write(join_paths(gname_, trim(adjustl(ci))), "oscstrr", oscstrr(:, o1, o2), [1], shape(oscstrr(:, o1, o2)))
        end do
      end do

      call h5%finalize()
    end subroutine write_excitons_hdf5
  
   
    subroutine hdf5_bandstructure_output(mpi_env, h5file, h5group, energies, energy_range, distances, label_names, label_distances, label_coordinates, characters)
      use xhdf5, only: xhdf5_type
      use modmpi, only: mpiinfo
      use precision, only: sp, dp
      use os_utils, only: join_paths

      type(mpiinfo), intent(in) :: mpi_env 
      character(*), intent(in) :: h5file, h5group

      real(dp), intent(in) :: energies(:, :), energy_range(2), distances(:), label_distances(:), label_coordinates(:, :)
      character(*), intent(in) :: label_names
      real(sp), intent(in), optional :: characters(:, :, :, :)

      ! Name constants
      character(*), parameter :: bands_group = 'bandstructure'
      character(*), parameter :: energies_dset = 'energies'
      character(*), parameter :: energy_range_dset  = 'energy_range'
      character(*), parameter :: distances_dset = 'distances'
      character(*), parameter :: label_names_dset = 'label_names'
      character(*), parameter :: label_distances_dset = 'label_distances'
      character(*), parameter :: label_coordinates_dset = 'label_coordinates'
      character(*), parameter :: characters_dset = 'characters'

      type(xhdf5_type) :: h5
      character(:), allocatable :: group

      call h5%initialize(h5file, mpi_env)
      call h5%initialize_group(h5group, bands_group)
      group = join_paths(h5group, bands_group)
      call h5%write(group, energies_dset, energies, [1, 1], shape(energies))
      call h5%write(group, energy_range_dset, energy_range, [1], shape(energy_range))
      call h5%write(group, distances_dset, distances, [1], shape(distances))
      call h5%write(group, label_names_dset, label_names)
      call h5%write(group, label_distances_dset, label_distances, [1], shape(label_distances))
      call h5%write(group, label_coordinates_dset, label_coordinates, [1, 1], shape(label_coordinates))
      if(present(characters)) then
        call h5%write(group, characters_dset, characters, [1, 1, 1, 1], shape(characters))
      end if
      call h5%finalize()
      
    end subroutine hdf5_bandstructure_output
end module m_write_hdf5
