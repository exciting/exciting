
subroutine plot_selfenergy()

    use modinput
    use modmain
    use modgw
    use mod_frequency
    use mod_hdf5
    use mod_mpi_gw
    use m_getunit
    use mod_pade
    implicit none
    ! local variables
    integer :: ik, ik_, ie, ie_, fid, recl
    integer :: npar, iom, nw
    character(20) :: s1, s2
    real(8) :: egap, enk, scr, sci, t1, df
    complex(8) :: ein, dsig, sigma, zt1
    
    integer :: nom, iom0(1), i
    real(8) :: wmin, wmax, div, eta
    real(8), allocatable :: om(:), w(:)
    
    ! parameters of the functions fitting the selfenergy
    complex(8), allocatable :: zx(:), zy(:)
    complex(8), allocatable :: a(:), sc(:), scw(:), poles(:)
    complex(8), allocatable :: selfc(:,:,:)
    real(8),    allocatable :: spectr(:,:,:)
    
    complex(8), allocatable :: p(:,:)

    complex(8), allocatable :: zvec(:)
    real(8),    allocatable :: tvec(:)
    
    integer :: n
    character(30) :: frmt1, frmt2, axis
       
    call init0
    input%groundstate%stypenumber = -1
    call init1

    nvelgw = chgval-occmax*dble(ibgw-1)
    nbandsgw = nbgw-ibgw+1
    call init_kqpoint_set
    call generate_freqgrid(freq, &
    &                      input%gw%freqgrid%fgrid, &
    &                      input%gw%freqgrid%fconv, &
    &                      input%gw%freqgrid%nomeg, &
    &                      input%gw%freqgrid%freqmax)
 
    if (allocated(evalsv)) deallocate(evalsv)
    allocate(evalsv(nstsv,kset%nkpt))
    allocate(vxcnn(ibgw:nbgw,kset%nkpt))
    call init_selfenergy(ibgw,nbgw,kset%nkpt,freq%nomeg)
    
    ! read data from files
#ifdef _HDF5_
    !call hdf5_read(fgwh5,"/","efermi",efermi)
    do ik = 1, kset%nkpt
      write(cik,'(I4.4)') ik
      path = "/kpoints/"//trim(adjustl(cik))
      call hdf5_read(fgwh5,path,"evalks",evalks(ibgw,ik),(/nbandsgw/))
      evalsv(:,ik) = evalks(:,ik)
      call hdf5_read(fgwh5,path,"vxcnn",vxcnn(ibgw,ik),(/nbandsgw/))
      call hdf5_read(fgwh5,path,"selfex",selfex(ibgw,ik),(/nbandsgw/))
      call hdf5_read(fgwh5,path,"selfec",selfec(ibgw,1,ik),(/nbandsgw,freq%nomeg/))
    end do ! ik
#else
    do ik = 1, kset%nkpt
      ik_ = kset%ikp2ik(ik)
      call getevalsvgw_new('GW_EVALSV.OUT',ik_,kqset%vkl(:,ik_), &
      &                     nstsv,evalsv(1,ik))
      evalks(ibgw:nbgw,ik) = evalsv(ibgw:nbgw,ik)
    end do
    call getunit(fid)
    open(fid,file='VXCNN.OUT',form='FORMATTED',status='OLD',action='READ')
    do ik = 1, kset%nkpt
      read(fid,*) s1, ik_, s2, kset%vkl(:,ik)
      do ie = ibgw, nbgw
        read(fid,*) ie_, vxcnn(ie,ik)
      end do
      read(fid,*)
    end do
    close(fid)
    if (allocated(evalqp)) deallocate(evalqp)
    allocate(evalqp(nstsv,kset%nkpt))
    call getevalqp(kset%nkpt,kset%vkl,evalqp)
    call readselfx
    call readselfc
#endif

    ! KS states analysis
    call fermi_exciting(input%groundstate%tevecsv, &
    &                   nvelgw, &
    &                   nbandsgw,kset%nkpt,evalks(ibgw:nbgw,:), &
    &                   kset%ntet,kset%tnodes,kset%wtet,kset%tvol, &
    &                   efermi,egap,df)
    call bandstructure_analysis('KS',ibgw,nbgw,kset%nkpt, &
    &                           evalks(ibgw:nbgw,:),efermi)

    ! QP states analysis
    call fermi_exciting(input%groundstate%tevecsv, &
    &                   nvelgw, &
    &                   nbandsgw,kset%nkpt,evalqp(ibgw:nbgw,:), &
    &                   kset%ntet,kset%tnodes,kset%wtet,kset%tvol, &
    &                   eferqp,egap,df)
    call bandstructure_analysis('G0W0',ibgw,nbgw,kset%nkpt, &
    &                           evalqp(ibgw:nbgw,:),eferqp)

    ! Frequency grid for the spectral function
    if (.not.associated(input%gw%selfenergy%SpectralFunctionPlot)) &
    & input%gw%selfenergy%SpectralFunctionPlot => getstructspectralfunctionplot(emptynode)
    nom = input%gw%selfenergy%SpectralFunctionPlot%nwgrid
    wmin = input%gw%selfenergy%SpectralFunctionPlot%wmin
    wmax = input%gw%selfenergy%SpectralFunctionPlot%wmax
    axis = trim(input%gw%selfenergy%SpectralFunctionPlot%axis)
    eta = input%gw%selfenergy%SpectralFunctionPlot%eta

    div = (wmax-wmin)/dble(nom)
    allocate(om(nom))
    om(:) = 0.d0  
    do iom = 1, nom
      om(iom) = wmin+dble(iom)*div
    end do
    iom0 = 0
    iom0 = minloc(om,om>0)
    ! write(*,*) 'iom0=', iom0
    ! write(*,*) om(iom0-1), om(iom0)

    ! doubled grid from -wmax to wmax
    allocate(zx(freq%nomeg))
    allocate(zy(freq%nomeg))
    do iom = 1, freq%nomeg
      zx(iom) = zi*freq%freqs(iom)
    end do

    ! Pade scheme only parameters 
    n = input%gw%selfenergy%npol
    allocate(p(n,n))

    ! MPF scheme only parameters 
    ! input%gw%selfenergy%npol = 2
    npar = 2*input%gw%selfenergy%npol
    allocate(a(npar))
    allocate(poles(npar))
    
    allocate(selfc(ibgw:nbgw,nom,kset%nkpt))
    selfc(:,:,:) = 0.d0
    allocate(spectr(ibgw:nbgw,nom,kset%nkpt))
    spectr(:,:,:) = 0.d0

    do ik = 1, kset%nkpt
      do ie = ibgw, nbgw

        ! enk = evalks(ie,ik)-efermi
        ! enk = evalqp(ie,ik)-eferqp
        do iom = 1, freq%nomeg
          zy(iom) = selfec(ie,iom,ik)
        end do

        if (iopac==0) then
          !-------------------
          ! Pade-approximant
          !-------------------
          select case(trim(axis))

            case('real')
              !----------------
              ! real frequency
              !----------------
              call padecof(freq%nomeg,zx,zy,n,p)
              do iom = 1, nom
                ein = cmplx(om(iom),eta,8)
                call gpade(freq%nomeg,zx,n,p,ein,sigma)
                if (om(iom)<0.d0) then
                  selfc(ie,iom,ik) = conjg(sigma)
                else
                  selfc(ie,iom,ik) = sigma
                end if
              end do
                        
            case('imag')
              !----------------
              ! complex frequency
              !----------------
              call padecof(freq%nomeg,zx,zy,n,p)
              do iom = 1, nom
                ein = cmplx(0.d0,om(iom),8)
                if (om(iom)<0.d0) then
                  call gpade(freq%nomeg,-zx,n,conjg(p),ein,sigma)
                else
                  call gpade(freq%nomeg,zx,n,p,ein,sigma)
                end if
                selfc(ie,iom,ik) = sigma
              end do

            case default
              stop "Error(plot_selfenergy): Unknown axis type!"

          end select

        else if (iopac==1) then
          !-------------------
          ! RGN Multipole-Fitting
          !-------------------
          select case(trim(axis))

            case('real')
              !----------------
              ! real frequency
              !----------------
              call setsac(iopac,freq%nomeg,npar,1.d0,selfec(ie,:,ik),freq%freqs,a,poles)
              do iom = 1, nom
                ein = cmplx(om(iom),eta,8)
                call getsac(iopac,nom,npar,1.d0,ein,om,a,sigma,dsig)
                if (om(iom)<0.d0) then
                  selfc(ie,iom,ik) = conjg(sigma)
                else
                  selfc(ie,iom,ik) = sigma
                end if
              end do
                    
            case('imag')
              !----------------
              ! complex frequency
              !----------------
              call setsac(iopac,freq%nomeg,npar,-1.d0,selfec(ie,:,ik),freq%freqs,a,poles)
              do iom = 1, iom0(1)-1
                ein = cmplx(0.d0,om(iom),8)
                call getsac(iopac,nom,npar,-1.d0,ein,om,a,sigma,dsig)
                selfc(ie,iom,ik) = sigma
              end do
              call setsac(iopac,freq%nomeg,npar,1.d0,selfec(ie,:,ik),freq%freqs,a,poles)
              do iom = iom0(1), nom
                ein = cmplx(0.d0,om(iom),8)
                call getsac(iopac,nom,npar,1.d0,ein,om,a,sigma,dsig)
                selfc(ie,iom,ik) = sigma
              end do

            case default
              stop "Error(plot_selfenergy): Unknown axis type!"

          end select

        end if ! AC type

        ! Spectral function
        do iom = 1, nom
          sigma = selfc(ie,iom,ik)
          scr = dble(selfex(ie,ik)+sigma)
          sci = aimag(selfex(ie,ik)+sigma)
          div = (om(iom)-enk-scr+dble(vxcnn(ie,ik)))**2+sci**2
          spectr(ie,iom,ik) = 1.d0/pi*dabs(sci)/div
        end do ! iom

      end do ! ie
    end do ! ik

    !---------
    ! ie = 4
    ! do iom = freq%nomeg, 1, -1
    !   write(7,'(3f18.6)') -freq%freqs(iom), conjg(selfec(ie,iom,1))
    ! end do
    ! do iom = 1, freq%nomeg
    !   write(7,'(3f18.6)') freq%freqs(iom), selfec(ie,iom,1)
    ! end do
    ! do iom = 1, nom
    !     write(8,'(3f18.6)') om(iom), selfc(ie,iom,1)
    ! end do
    !---------
    
    ik = 1
    kset%nkpt = 1

    n = nbgw-ibgw+1
    write(frmt1,'("(",i8,"f14.6)")') 1+2*n
    write(*,*) trim(frmt1)
    write(frmt2,'("(",i8,"f14.6)")') 1+n
    write(*,*) trim(frmt2)
    
    ! OUTPUT: integrated over BZ
    open(73,file='SelfC-Re.dat',form='FORMATTED',status='UNKNOWN',action='WRITE')
    open(74,file='SelfC-Im.dat',form='FORMATTED',status='UNKNOWN',action='WRITE')
    open(75,file='SelfXC.dat',form='FORMATTED',status='UNKNOWN',action='WRITE')
    open(76,file='SpectralFunction.dat',form='FORMATTED',status='UNKNOWN',action='WRITE')
    open(77,file='GF.dat',form='FORMATTED',status='UNKNOWN',action='WRITE')
    allocate(tvec(ibgw:nbgw),zvec(ibgw:nbgw))
    do ik = 1, kset%nkpt
      write(74,*) '# ik = ', ik
      write(75,*) '# ik = ', ik
      write(76,*) '# ik = ', ik
      write(77,*) '# ik = ', ik
      do iom = 1, nom
        do ie = ibgw, nbgw
          zvec(ie) = selfex(ie,ik)+selfc(ie,iom,ik)
          tvec(ie) = om(iom)-enk-dble(zvec(ie))+dble(vxcnn(ie,ik))
        end do
        ! \Sigma_c
        write(73,trim(frmt1)) om(iom), dble(selfc(:,iom,ik))
        write(74,trim(frmt1)) om(iom), imag(selfc(:,iom,ik))
        ! \Sigma_xc
        write(75,trim(frmt1)) om(iom), zvec
        ! Spectral function
        write(76,trim(frmt2)) om(iom), spectr(:,iom,ik)
        ! Denominator
        write(77,trim(frmt2)) om(iom), tvec 
      end do
      write(74,*)
      write(75,*)
      write(76,*)
      write(77,*)
    end do
    deallocate(tvec,zvec)
    close(74)
    close(75)
    close(76)
    close(77)
    
    ! clear memory
    deallocate(om)
    deallocate(a,poles)
    deallocate(p)
    call delete_selfenergy
    deallocate(evalsv)
    deallocate(vxcnn)
    deallocate(spectr)
    deallocate(selfc)
      
    return
end subroutine

