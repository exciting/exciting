! Copyright(C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_write_hdf5

  implicit none
  ! filename for intermediate HDF5 output, i.e. the BSE matrix elements.
  character(256) :: fhdf5_inter
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
      character(256) :: gname_, group, momentum_index, epsname, ci
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
      write(momentum_index, '(I4.4)') iq ! Generate string out of momentum transfer index
      if (.not. hdf5_exist_group(fhdf5, gname_, momentum_index )) then
        call hdf5_create_group(fhdf5,gname_, momentum_index)
      end if
      ! Create parameter Subgroup
      group=trim(gname_)//trim(momentum_index)
      if (.not. hdf5_exist_group(fhdf5, group, "parameters")) then
        call hdf5_create_group(fhdf5,group, "parameters")
      end if
      ! Write meta data
      group=trim(gname_)//trim(adjustl(momentum_index))//"/parameters"     
      call hdf5_write(fhdf5,group,"ivgmt", ivgmt(1,iq), shape(ivgmt(1:3,iq)))
      call hdf5_write(fhdf5,group,"vqlmt",vqlmt(1,iq), shape(vqlmt(1:3,iq)))
      call hdf5_write(fhdf5,group,"vgcmt",vgcmt(1,iq), shape(vgcmt(1:3,iq)))
      call hdf5_write(fhdf5,group,"vqcmt",vqcmt(1,iq), shape(vqcmt(1:3,iq)))
      call hdf5_write(fhdf5,group,"escale",escale)
      call hdf5_write(fhdf5,group,"broad",escale*input%xs%broad)
      call hdf5_write(fhdf5,group,"nk_bse",nk_bse)
      call hdf5_write(fhdf5,group,"foff",foff)
      
      ! Write dielectric function
      group=trim(gname_)//trim(adjustl(momentum_index))
      if (.not. hdf5_exist_group(fhdf5, group, "diel" )) then
        call hdf5_create_group(fhdf5,group, "diel")
      end if
      group=trim(gname_)//trim(adjustl(momentum_index))//"/diel"
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
      group=trim(gname_)//trim(adjustl(momentum_index))
      if (.not. hdf5_exist_group(fhdf5, group, "loss" )) then
        call hdf5_create_group(fhdf5,group, "loss")
      end if
      group=trim(gname_)//trim(adjustl(momentum_index))//"/loss"
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
      group=trim(gname_)//trim(adjustl(momentum_index))
      if (.not. hdf5_exist_group(fhdf5, group, "sigma" )) then
        call hdf5_create_group(fhdf5,group, "sigma")
      end if
      group=trim(gname_)//trim(adjustl(momentum_index))//"/sigma"
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
      use modxs, only: escale, unitout, unit1
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
      real(8) :: pm
      character(256) :: fnexc, frmt, tdastring, bsetypestring, tistring, scrtypestring
      character(256) :: syscommand, excitondir
      character(128) :: gname_, group, ci, momentum_index
      
#ifdef _HDF5_
      ! Create group excitons-bsetypestring-scrtypestring
      if (.not. hdf5_exist_group(fhdf5, "/", gname)) then
        call hdf5_create_group(fhdf5,"/", gname)
      end if
      gname_="/"//trim(adjustl(gname))//"/"
      
      ! If necessary create subgroup for each momentum transfer index
      if (present(iqmt)) then
        write(momentum_index, '(I4.4)') iqmt
        if (.not. hdf5_exist_group(fhdf5, gname_, momentum_index)) then
          call hdf5_create_group(fhdf5,gname_, momentum_index)
        end if
        gname_="/"//trim(adjustl(gname))//"/"//trim(adjustl(momentum_index))//"/"
      end if
        
      
      ! Write Meta Data
      if (.not. hdf5_exist_group(fhdf5, gname_, "parameters")) then
        call hdf5_create_group(fhdf5,gname_, "parameters")
      end if
      iq=iqmt
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
        evalre_(lambda)=evalre(lambda)*escale
      end do
      call hdf5_write(fhdf5,gname_,"evalre", evalre_(1), shape(evalre_))
      deallocate(evalre_)
      

      io1=1
      io2=1
      if(iq == 1) io2 = 3
      
      do o1=io1,io2
        write(ci, '(I4.1)') o1
        gname_="/"//trim(adjustl(gname))//"/"//trim(adjustl(momentum_index))//"/"
        if (.not. hdf5_exist_group(fhdf5, gname_, trim(adjustl(ci)))) then
          call hdf5_create_group(fhdf5,gname_,trim(adjustl(ci)))
        end if
        group=trim(adjustl(gname_))//trim(adjustl(ci))//"/"
        call hdf5_write(fhdf5,group,"oscstrr", oscstrr(1,o1), shape(oscstrr(1:nexc,o1)))
      end do
#endif
    end subroutine write_excitons_hdf5
    
    subroutine write_weights_hdf5(lambda ,vkl, vkl0, ivmin, ivmax, icmin, icmax, rv, rc, arv, arc)
      use modinput, only: input  
      use mod_hdf5

      implicit none
      ! excitonic index
      integer, intent(in) :: lambda
      ! k-vectors in lattice coordinates
      real(8), intent(in) :: vkl(:,:), vkl0(:,:)
      ! min and max values for the valence and conduction bands 
      integer, intent(in) ::  ivmin, ivmax, icmin, icmax
      ! resonant valence and conduction band weights
      real(8), intent(in) :: rv(:,:), rc(:,:)
      ! for non-TDA calculations, the anti-resonant weights are stored as well
      real(8), intent(in), optional :: arv(:,:), arc(:,:)
      ! local variables
      character(256) :: gname_, bsetypestring, scrtypestring, tdastring, pos_, lambda_, params_
      real(8), allocatable :: inter(:,:)
#ifdef _HDF5_
      ! determine the name of the group
      ! determine the TDA string
      if (input%xs%bse%coupling) then
        tdastring=''
      else
        if (input%xs%bse%chibarq) then
          tdastring='-TDA-BAR'
        else
          tdastring='-TDA'
        end if
      end if
      ! determine bsetypestring & scrtypestring
      bsetypestring = '-'//trim(input%xs%bse%bsetype)//trim(tdastring)
      scrtypestring = '-'//trim(input%xs%screening%screentype)
      gname_='weights'//trim(bsetypestring)//trim(scrtypestring)
      ! generate group
      if (.not. hdf5_exist_group(fhdf5,'/', gname_)) then
        call hdf5_create_group(fhdf5,'/', gname_)
      end if
      gname_='/'//trim(gname_)//'/'
      if (.not. hdf5_exist_group(fhdf5,gname_, 'parameters')) then
        params_=trim(adjustl(gname_))//'parameters/'
        call hdf5_create_group(fhdf5,gname_, 'parameters')
        call hdf5_write(fhdf5,params_,"vkl", vkl(1,1), shape(vkl(:,:)))
        call hdf5_write(fhdf5,params_,"vkl0", vkl0(1,1), shape(vkl0(:,:)))
        call hdf5_write(fhdf5,params_,"ivmin",ivmin)
        call hdf5_write(fhdf5,params_,"ivmax",ivmax)
        call hdf5_write(fhdf5,params_,"icmin",ivmin)
        call hdf5_write(fhdf5,params_,"icmax",ivmin)
      end if
      if (.not. hdf5_exist_group(fhdf5,gname_, 'rvwgrid')) then
        call hdf5_create_group(fhdf5,gname_, 'rvwgrid')
      end if
      if (.not. hdf5_exist_group(fhdf5,gname_, 'rcwgrid')) then
        call hdf5_create_group(fhdf5,gname_, 'rcwgrid')
      end if
      if ((present(arv)) .and. (.not. hdf5_exist_group(fhdf5,gname_, 'arvwgrid'))) then
        call hdf5_create_group(fhdf5,gname_, 'arvwgrid')
      end if
      if ((present(arc)) .and. (.not. hdf5_exist_group(fhdf5,gname_, 'arvcgrid'))) then
        call hdf5_create_group(fhdf5,gname_, 'arcwgrid')
      end if
      write(lambda_,'(I4.4)') lambda
      pos_=trim(adjustl(gname_))//'rvwgrid/'
      call hdf5_write(fhdf5,pos_,lambda_,rv(1,1), shape(rv(:,:)))
      pos_=trim(adjustl(gname_))//'rcwgrid/'
      call hdf5_write(fhdf5,pos_,lambda_,rc(1,1), shape(rc(:,:)))
      if (present(arv)) then
        pos_=trim(adjustl(gname_))//'arvwgrid/'
        call hdf5_write(fhdf5,pos_,lambda_,arv(1,1), shape(arv(:,:)))
      end if
      if (present(arc)) then
        pos_=trim(adjustl(gname_))//'arcwgrid/'
        call hdf5_write(fhdf5,pos_,lambda_,arc(1,1), shape(arc(:,:)))
      end if
#endif
    end subroutine write_weights_hdf5
    
    subroutine write_kpathplot_hdf5(lambda,iv1,iv2,ic1,ic2,rvw,rcw,arvw,arcw)
      use mod_hdf5
      use m_read_bandstructure, only: kpathlength_, energyval_
      use modxs, only: escale
      use modinput, only: input

      implicit none

      integer, intent(in) :: lambda, iv1, iv2, ic1, ic2
      real(8), intent(in) :: rvw(:,:), rcw(:,:)
      real(8), intent(in), optional :: arvw(:,:), arcw(:,:)
      ! local variables
      character(256) :: gname_, params_, pos_, bsetypestring, scrtypestring, &
        &               lambda_, tdastring

#ifdef _HDF5_
      ! create group for kpathweights 
      if (input%xs%bse%coupling) then
        tdastring=''
      else
        if (input%xs%bse%chibarq) then
          tdastring='-TDA-BAR'
        else
          tdastring='-TDA'
        end if
      end if
      bsetypestring = '-'//trim(input%xs%bse%bsetype)//trim(tdastring)
      scrtypestring = '-'//trim(input%xs%screening%screentype)
      gname_='kpathweights'//trim(bsetypestring)//trim(scrtypestring)
      if (.not. hdf5_exist_group(fhdf5,'/', gname_)) then
        call hdf5_create_group(fhdf5,'/', gname_)
      endif
      gname_='/'//trim(adjustl(gname_))//'/'
      
      ! create parameter group and write parameters
      if (.not. hdf5_exist_group(fhdf5,gname_, 'parameters')) then
        call hdf5_create_group(fhdf5,gname_, 'parameters')
        params_=trim(adjustl(gname_))//'parameters/'
        
        call hdf5_write(fhdf5,params_,"iv1", iv1)
        call hdf5_write(fhdf5,params_,"iv2", iv2)
        call hdf5_write(fhdf5,params_,"ic1", ic1)
        call hdf5_write(fhdf5,params_,"ic2", ic2)
        call hdf5_write(fhdf5,params_,"escale", escale)
        call hdf5_write(fhdf5, params_, 'kpathlength',kpathlength_(1,1), shape(kpathlength_))
        call hdf5_write(fhdf5, params_, 'energyval_',energyval_(1,1), shape(energyval_))
      endif

      ! create groups for valence and conduction weights
      if (.not. hdf5_exist_group(fhdf5,gname_, 'rvw')) then
        call hdf5_create_group(fhdf5,gname_, 'rvw')
      end if
      if (.not. hdf5_exist_group(fhdf5,gname_, 'rcw')) then
        call hdf5_create_group(fhdf5,gname_, 'rcw')
      end if
      if ((present(arvw)) .and. (.not. hdf5_exist_group(fhdf5,gname_, 'arvw'))) then
        call hdf5_create_group(fhdf5,gname_, 'arvw')
      end if
      if ((present(arcw)) .and. (.not. hdf5_exist_group(fhdf5,gname_, 'arcw'))) then
        call hdf5_create_group(fhdf5,gname_, 'arvw')
      end if
      ! write actual data
      write(lambda_,'(I4.4)') lambda
      pos_=trim(adjustl(gname_))//'rvw'
      call hdf5_write(fhdf5,pos_,lambda_,rvw(1,1), shape(rvw(:,:)))
      pos_=trim(adjustl(gname_))//'rcw'
      call hdf5_write(fhdf5,pos_,lambda_,rcw(1,1), shape(rcw(:,:)))
      if (present(arvw)) then
        pos_=trim(adjustl(gname_))//'arvw'
        call hdf5_write(fhdf5,pos_,lambda_,arvw(1,1), shape(arvw(:,:)))
      endif
      if (present(arcw)) then
        pos_=trim(adjustl(gname_))//'arcw'
        call hdf5_write(fhdf5,pos_,lambda_,arcw(1,1), shape(arcw(:,:)))
      endif
#endif
    end subroutine write_kpathplot_hdf5

    subroutine write_bandstr_hdf5(nkpt,nstsv,dpp1d,dvp1d,evalsv,bc)
      use mod_hdf5
      use mod_kpoint, only: vkl
      use modmain, only: natoms 
      use mod_atoms, only: idxas, nspecies
      use modinput, only: input

      implicit none
      integer, intent(in) :: nkpt, nstsv
      real(8), intent(in) :: evalsv(:,:)
      real(8), intent(in) :: dpp1d(:)
      real(8), intent(in) :: dvp1d(:)
      real(4), intent(in), optional :: bc(:,:,:,:)
      ! local variables 
      character(256) :: gname_, params_
      character(256), allocatable :: vertices_(:)
      integer :: shape_(1), iv
      ! create array of labels
      if (allocated(vertices_)) deallocate(vertices_)
      shape_=shape(dvp1d)
      allocate(vertices_(shape_(1)))
      do iv=1, shape_(1)
        vertices_(iv)=trim(input%properties%bandstructure%plot1d%path%pointarray(iv)%point%label)
      end do
      ! create bandstructure group
      if (.not. hdf5_exist_group(fhdf5,'/', 'bandstructure')) then
        call hdf5_create_group(fhdf5,'/', 'bandstructure')
      end if 
      gname_='/bandstructure/'
      if (.not. hdf5_exist_group(fhdf5,gname_, 'parameters')) then
        call hdf5_create_group(fhdf5, gname_, 'parameters')
        params_=trim(adjustl(gname_))//'parameters/'
        if (present(bc)) then
          call hdf5_write(fhdf5,params_,'idxas',idxas(1,1),shape(idxas))
          call hdf5_write(fhdf5,params_,'nspecies',nspecies)
          call hdf5_write(fhdf5,params_,'natoms',natoms(1),shape(natoms))
        end if
        call hdf5_write(fhdf5,params_,'vkl',vkl(1,1),shape(vkl))
        call hdf5_write(fhdf5,params_,'nkpt',nkpt)
        call hdf5_write(fhdf5,params_,'nstsv',nstsv)
        call hdf5_write(fhdf5, params_, 'labelpoints', dvp1d(1), shape(dvp1d))
        call hdf5_write(fhdf5, params_, 'vertex', vertices_(1), shape(vertices_))
      end if
      ! write data
      call hdf5_write(fhdf5, gname_, 'points', dpp1d(1), shape(dpp1d))
      call hdf5_write(fhdf5, gname_, 'evalsv', evalsv(1,1), shape(evalsv))
      if (present(bc)) then
        call hdf5_write(fhdf5, gname_, 'bc', bc(1,1,1,1), shape(bc))
      end if
#ifdef _HDF5_
      
#endif
    end subroutine write_bandstr_hdf5
end module m_write_hdf5
