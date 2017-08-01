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
      & gname, evalim, oscstra, sort, iqmt)
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
      real(8), intent(in), optional :: evalim(hamsize)
      complex(8), intent(in), optional :: oscstra(:,:)
      logical, intent(in), optional :: sort
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
      
      if(present(sort)) then 
        fsort = sort
      else
        fsort = .false.
      end if
      if (present(iqmt)) then
        iq=iqmt
      else
        iq=1
      end if

#ifdef _HDF5_
      allocate(idxsort(hamsize))
      if(fsort) then 
        allocate(idxsort2(hamsize/2))
        allocate(evalre_sorted(hamsize))
        call sortidx(hamsize, evalre, idxsort)
        evalre_sorted = evalre(idxsort)
        idxsort = cshift(idxsort, hamsize/2)
        evalre_sorted = evalre(idxsort)
        do i = 1, hamsize/2
          idxsort2(i) = idxsort(hamsize-i+1)
        end do
        idxsort(hamsize/2+1:hamsize) = idxsort2
        evalre_sorted = evalre(idxsort)
        deallocate(idxsort2, evalre_sorted)
      else
        do lambda = 1, nexc
          idxsort(lambda) = lambda
        end do
      end if
      
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
      group=trim(adjustl(gname_))//"parameters"     
      call hdf5_write(fhdf5,group,"ivgmt", ivgmt(1,iq), shape(ivgmt(1:3,iq)))
      call hdf5_write(fhdf5,group,"vqlmt",vqlmt(1,iq), shape(vqlmt(1:3,iq)))
      call hdf5_write(fhdf5,group,"vgcmt",vgcmt(1,iq), shape(vgcmt(1:3,iq)))
      call hdf5_write(fhdf5,group,"vqcmt",vqcmt(1,iq), shape(vqcmt(1:3,iq)))
      call hdf5_write(fhdf5,group,"escale",escale)
      call hdf5_write(fhdf5,group,"eshift",eshift*escale)
      
      ! Write real part of energy eigenvalues
      allocate(evalre_(nexc))
      do lambda=1, nexc
        evalre_(lambda)=evalre(idxsort(lambda))*escale
      end do
      call hdf5_write(fhdf5,gname_,"evalre", evalre_(1), shape(evalre_))
      deallocate(evalre_)
      
      if (present(evalim)) then
        allocate(evalim_(nexc))
        do lambda=1,nexc
          evalim_(lambda)=evalim(idxsort(lambda))*escale
        end do
        call hdf5_write(fhdf5,gname_,"evalim", evalim_(1), shape(evalim_))
        deallocate(evalim_)
      end if

      ! Loop over optical components
      if(present(iqmt)) then 
        iq = iqmt
      else
        iq = 1
      end if

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
        allocate(oscstrr_(nexc))
        do lambda=1, nexc
          oscstrr_(lambda)=oscstrr(idxsort(lambda),o1)
        end do
        call hdf5_write(fhdf5,group,"oscstrr", oscstrr_(1), shape(oscstrr_))
        deallocate(oscstrr_)

        if (present(oscstra)) then
        allocate(oscstra_(nexc))
        do lambda=1, nexc
          oscstra_(lambda)=oscstra(idxsort(lambda),o1)
        end do
        call hdf5_write(fhdf5,group,"oscstra", oscstra_(1), shape(oscstra_))
        deallocate(oscstra_)
        end if

      end do
#endif
    end subroutine write_excitons_hdf5

end module m_write_hdf5
