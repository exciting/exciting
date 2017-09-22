! Copyright(C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_write_hdf5

  implicit none

  contains

    subroutine write_spectra_hdf5(iq, foff, w, eps, loss, sigma, gname)
      use mod_lattice, only: omega
      use mod_misc, only: version
      use modinput
      use modmpi
      use modbse, only: nk_bse
      use modxs, only: escale, unitout, unit1, fhdf5
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
      character(256) :: gname_, group, ciq, epsname, ci
      complex(8), allocatable :: eps_(:)
      real(8), allocatable :: loss_(:)


      !Call kramkron(iop1, iop2, 1.d-8, n, w, imeps, kkeps)

      igqmt = ivgigq(ivgmt(1,iq),ivgmt(2,iq),ivgmt(3,iq),iq)
#ifdef _HDF5_
      ! Create Group 'spectra-bsetypestring-scrtypestring'
      if (.not. hdf5_exist_group(fhdf5, "/", gname)) then
        call hdf5_create_group(fhdf5,"/", gname)
      end if
      gname_="/"//trim(adjustl(gname))//"/"
      ! Create Subgroup for each momentum transfer entry
      write(ciq, '(I4.4)') iq ! Generate string out of momentum transfer index
      if (.not. hdf5_exist_group(fhdf5, gname_, ciq )) then
        call hdf5_create_group(fhdf5,gname_, ciq)
      end if
      ! Create parameter Subgroup
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
      if (foff) then ! write full dielectric tensor
        call hdf5_write(fhdf5,group,"epsm", eps(1,1,1), shape(eps(:,:,:)))
      else ! write the diagonal entries of the dielectric tensor separately
        do i=1,3
          write(ci, '(I4.2)') i*10+i
          epsname='epsm('//trim(adjustl(ci))//')'
          allocate(eps_(size(w)))
          eps_(:)=eps(i,i,:)
          call hdf5_write(fhdf5,group,epsname, eps_(1), shape(eps_(:)))
          deallocate(eps_)
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
          allocate(loss_(size(w)))
          loss_(:)=loss(i,i,:)
          epsname='lossfct('//trim(adjustl(ci))//')' 
          call hdf5_write(fhdf5,group,epsname, loss_(1), shape(loss_(:)))
          deallocate(loss_)
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

#endif   
   end subroutine write_spectra_hdf5

   subroutine write_excitons_hdf5(hamsize, nexc, eshift, evalre, oscstrr,&
      & gname, iqmt)
      use mod_lattice, only: omega
      use mod_misc, only: version
      use modinput
      use modmpi
      use modbse, only: nk_bse
      use modxs, only: escale, unitout, unit1, fhdf5
      use modxs, only: ivgmt, vqlmt, vgcmt, vqcmt
      use modxs, only: sptclg, ivgigq
      use fox_wxml
      use m_getunit
      use mod_hdf5

      implicit none

      ! Arguments
      integer(4), intent(in) :: hamsize, nexc
      real(8), intent(in) :: eshift
      real(8), intent(in) :: evalre(hamsize)
      complex(8), intent(in) :: oscstrr(:,:)
      character(128), intent(in) :: gname
      integer(4), intent(in), optional :: iqmt
   
      ! Local
      logical :: fsort
      integer(4) :: o1, lambda, unexc, i, io1, io2, iq
      integer(4), allocatable :: idxsort(:), idxsort2(:)
      real(8), allocatable :: evalre_sorted(:)
      real(8), allocatable :: evalim_(:), evalre_(:)
      complex(8), allocatable :: oscstrr_(:), oscstra_(:)
      real(8) :: pm
      character(256) :: fnexc, frmt, tdastring, bsetypestring, tistring, scrtypestring
      character(256) :: syscommand, excitondir
      character(128) :: gname_, group, ci, ciq
#ifdef DGRID
      character(256) :: dgrid_dotext
#endif
      
#ifdef _HDF5_
      ! Create group excitons-bsetypestring-scrtypestring
      if (.not. hdf5_exist_group(fhdf5, "/", gname)) then
        call hdf5_create_group(fhdf5,"/", gname)
      end if
      gname_="/"//trim(adjustl(gname))//"/"
      
      ! If necessary create subgroup for each momentum transfer index
      if (present(iqmt)) then
        write(ciq, '(I4.4)') iqmt
        if (.not. hdf5_exist_group(fhdf5, gname_, ciq)) then
          call hdf5_create_group(fhdf5,gname_, ciq)
        end if
        gname_="/"//trim(adjustl(gname))//"/"//trim(adjustl(ciq))//"/"
      end if
        
      
      ! Write Meta Data
      if (.not. hdf5_exist_group(fhdf5, gname_, "parameters")) then
        call hdf5_create_group(fhdf5,gname_, "parameters")
      end if
      print *, 'iqmt=', iqmt
      iq=iqmt
      group=trim(adjustl(gname_))//"parameters"
      print *, 'group=', group     
      call hdf5_write(fhdf5,group,"ivgmt", ivgmt(1,iq), shape(ivgmt(1:3,iq)))
      call hdf5_write(fhdf5,group,"vqlmt",vqlmt(1,iq), shape(vqlmt(1:3,iq)))
      call hdf5_write(fhdf5,group,"vgcmt",vgcmt(1,iq), shape(vgcmt(1:3,iq)))
      call hdf5_write(fhdf5,group,"vqcmt",vqcmt(1,iq), shape(vqcmt(1:3,iq)))
      call hdf5_write(fhdf5,group,"escale",escale)
      call hdf5_write(fhdf5,group,"eshift",eshift*escale)
      
      ! Write real part of energy eigenvalues
      allocate(evalre_(nexc))
      do lambda=1, nexc
        evalre_(lambda)=evalre(lambda)*escale
      end do
      call hdf5_write(fhdf5,gname_,"evalre", evalre_(1), shape(evalre_))
      deallocate(evalre_)
      

      io1=1
      io2=3
      if(iq /= 1) then 
        io1=1
        io2=1
      end if
      
      do o1=io1,io2
        write(ci, '(I4.2)') o1*10+o1
        write(*,*) 'ci in write_exciton:', ci
        gname_="/"//trim(adjustl(gname))//"/"//trim(adjustl(ciq))//"/"
        if (.not. hdf5_exist_group(fhdf5, gname_, trim(adjustl(ci)))) then
          call hdf5_create_group(fhdf5,gname_,trim(adjustl(ci)))
        end if
        group=trim(adjustl(gname_))//trim(adjustl(ci))//"/"
        write(*,*) 'group in write_exciton:', group
        call hdf5_write(fhdf5,group,"oscstrr", oscstrr(1,o1), shape(oscstrr(1:nexc,o1)))
      end do
#endif
    end subroutine write_excitons_hdf5

end module m_write_hdf5
