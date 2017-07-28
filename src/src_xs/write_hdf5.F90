! Copyright(C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_write_hdf5

  implicit none

  contains

    subroutine write_hdf5(iq, foff, w, eps, loss, sigma, gname)
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
      character(256) :: fhdf5, gname_, group, ciq, epsname, ci


      !Call kramkron(iop1, iop2, 1.d-8, n, w, imeps, kkeps)

      igqmt = ivgigq(ivgmt(1,iq),ivgmt(2,iq),ivgmt(3,iq),iq)
      fhdf5="bse_output.h5"
#ifdef _HDF5_
      call hdf5_initialize()
      call hdf5_create_file(fhdf5)
      if (.not. hdf5_exist_group(fhdf5, "/", gname)) then
        call hdf5_create_group(fhdf5,"/", gname)
      end if
      gname_="/"//trim(adjustl(gname))//"/"
      write(ciq, '(I4.4)') iq ! Generate string out of momentum transfer index
      if (.not. hdf5_exist_group(fhdf5, gname_, ciq )) then
        call hdf5_create_group(fhdf5,gname_, ciq)
      end if
      group=trim(gname_)//trim(ciq)
      if (.not. hdf5_exist_group(fhdf5, group, "parameters")) then
        call hdf5_create_group(fhdf5,group, "parameters")
      end if
      ! Write meta data
      group=trim(gname_)//trim(adjustl(ciq))//"/parameters"     
      call hdf5_write(fhdf5,group,"ivgmt", ivgmt(1,iq), shape(ivgmt(1:3,iq)))
      call hdf5_write(fhdf5,group,"vqlmt",vqlmt(1,iq), shape(vqlmt(1:3,iq)))
      call hdf5_write(fhdf5,group,"vgcmt",vgcmt(1,iq), shape(vgcmt(1:3,iq)))
      call hdf5_write(fhdf5,group,"vqcmt",vqcmt(1,iq), shape(vqcmt(1:3,iq)))
      call hdf5_write(fhdf5,group,"escale",escale)
      call hdf5_write(fhdf5,group,"broad",escale*input%xs%broad)
      call hdf5_write(fhdf5,group,"nk_bse",nk_bse)
      call hdf5_write(fhdf5,group,"foff",foff)
      
      ! Write dielectric function
      group=trim(gname_)//trim(adjustl(ciq))
      if (.not. hdf5_exist_group(fhdf5, group, "diel" )) then
        call hdf5_create_group(fhdf5,group, "diel")
      end if
      group=trim(gname_)//trim(adjustl(ciq))//"/diel"
      do iw=1, size(w)
          w_(iw)=w(iw)*escale
      end do
      call hdf5_write(fhdf5,group,"w", w_(1), shape(w_(:)))
      write(*,*) 'foff=', foff
      write(*,*) input%xs%dfoffdiag
      if (foff) then ! write full dielectric tensor
        call hdf5_write(fhdf5,group,"epsm", eps(1,1,1), shape(eps(:,:,:)))
      else ! write the diagonal entries of the dielectric tensor separately
        do i=1,3
          write(ci, '(I4.2)') i*10+i
          epsname='epsm('//trim(adjustl(ci))//')' 
          call hdf5_write(fhdf5,group,epsname, eps(i,i,1), shape(eps(i,i,:)))
        end do
      end if
     
      ! Write loss function
      group=trim(gname_)//trim(adjustl(ciq))
      if (.not. hdf5_exist_group(fhdf5, group, "loss" )) then
        call hdf5_create_group(fhdf5,group, "loss")
      end if
      group=trim(gname_)//trim(adjustl(ciq))//"/loss"
      call hdf5_write(fhdf5,group,"w", w_(1), shape(w_(:)))
      if (foff) then ! write full loss tensor
        call hdf5_write(fhdf5,group,"lossfct", loss(1,1,1), shape(loss(:,:,:)))
      else ! write the diagonal entries of the loss function separately
        do i=1,3
          write(ci, '(I4.2)') i*10+i
          epsname='lossfct('//trim(adjustl(ci))//')' 
          call hdf5_write(fhdf5,group,epsname, loss(i,i,1), shape(loss(i,i,:)))
        end do
      end if
      
      ! Write sigma
      group=trim(gname_)//trim(adjustl(ciq))
      if (.not. hdf5_exist_group(fhdf5, group, "sigma" )) then
        call hdf5_create_group(fhdf5,group, "sigma")
      end if
      group=trim(gname_)//trim(adjustl(ciq))//"/sigma"
      call hdf5_write(fhdf5,group,"w", w_(1), shape(w_(:)))
      call hdf5_write(fhdf5,group,"sigma", sigma(1), shape(sigma(:)))

      call hdf5_finalize()
#endif    
   end subroutine write_hdf5

end module m_write_hdf5
