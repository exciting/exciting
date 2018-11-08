
subroutine spectralFunction()
    use modinput
    use modmain
    use modgw
    use mod_frequency
    use mod_vxc
    use mod_aaa_approximant
    use mod_pade
    implicit none
    ! local variables
    type(aaa_approximant) :: aaa_minus, aaa_plus
    integer(4) :: nw, iw, iom
    real(8)    :: wmin, wmax, div, eta, tol
    real(8)    :: sRe, sIm, enk
    complex(8) :: sc, zw, dsc
    
    integer(4) :: ik, ie, n
    character(80) :: axis, frmt1, frmt2
        
    ! parameters of the functions fitting the selfenergy
    real(8),    allocatable :: w(:), sf(:,:,:)
    complex(8), allocatable :: zj(:), fj(:)
    complex(8), allocatable :: sc_ac(:,:,:)

    !-----------------
    ! Initialization
    !-----------------
    call init0
    input%groundstate%stypenumber = -1
    call init1
    call init_kqpoint_set
 
    nvelgw = chgval-occmax*dble(ibgw-1)
    nbandsgw = nbgw-ibgw+1

    allocate(vxcnn(ibgw:nbgw,kset%nkpt))
    call read_vxcnn()
    
    call init_selfenergy(ibgw,nbgw,kset%nkpt)
    call readselfx()
    call readselfc()
    call readevalqp()

    ! Frequency grid for the spectral function
    if ( .not.associated(input%gw%selfenergy%SpectralFunctionPlot) ) &
        input%gw%selfenergy%SpectralFunctionPlot => getstructspectralfunctionplot(emptynode)
    nw   = input%gw%selfenergy%SpectralFunctionPlot%nwgrid
    wmin = input%gw%selfenergy%SpectralFunctionPlot%wmin
    wmax = input%gw%selfenergy%SpectralFunctionPlot%wmax
    axis = trim(input%gw%selfenergy%SpectralFunctionPlot%axis)
    eta  = input%gw%selfenergy%SpectralFunctionPlot%eta
    tol  = input%gw%selfenergy%tol

    if (axis /= 'real' .and. axis /= 'imag') stop 'Error(spectralFunction): Wrong axis type!'
    
    allocate( w(nw) )
    w(:) = 0.d0  
    div  = (wmax-wmin) / dble(nw-1)
    do iw = 1, nw
        w(iw) = wmin + dble(iw-1)*div
    end do

    allocate(zj(freq_selfc%nomeg), fj(freq_selfc%nomeg))

    allocate(sc_ac(ibgw:nbgw,nw,kset%nkpt))
    sc_ac(:,:,:) = zzero
    allocate(sf(ibgw:nbgw,nw,kset%nkpt))
    sf(:,:,:) = zzero

    do ik = 1, kset%nkpt
        do ie = ibgw, nbgw

            do iom = 1, freq_selfc%nomeg
                zj(iom) = cmplx(0.d0, freq_selfc%freqs(iom), 8)
                fj(iom) = selfec(ie,iom,ik)
            end do
            
            if (input%gw%selfenergy%actype == 'aaa' ) then
                ! Compute the approximant
                if (axis == 'real') then
                    call set_aaa_approximant(aaa_plus, -zj, fj, tol)  ! <----- No idea why only this works ????
                else
                    call set_aaa_approximant(aaa_plus, zj, fj, tol)
                end if
                call init_aaa_approximant(aaa_minus, aaa_plus%nj, -aaa_plus%zj, &
                                          conjg(aaa_plus%fj), conjg(aaa_plus%wj))
            end if

            enk = evalks(ie,ik)-eferks

            do iw = 1, nw
            
                if (axis == 'real') then
                    zw = cmplx(w(iw), 0.d0, 8)
                else
                    zw = cmplx(0.d0, w(iw), 8)
                end if

                if (w(iw) < 0.d0) then

                    if (input%gw%selfenergy%actype == 'aaa') then
                        sc = get_aaa_approximant(aaa_minus, zw)
                    else if (input%gw%selfenergy%actype == 'pade') then
                        call pade_approximant(size(zj), -zj, conjg(fj), zw, sc, dsc)
                    end if
                    
                else
                    
                    if (input%gw%selfenergy%actype == 'aaa') then
                        sc = get_aaa_approximant(aaa_plus, zw)
                    else if (input%gw%selfenergy%actype == 'pade') then
                        call pade_approximant(size(zj), zj, fj, zw, sc, dsc)
                    end if

                end if

                ! AC correlation self-energy
                sc_ac(ie,iw,ik) = sc

                ! Spectral function
                dsc = selfex(ie,ik) + sc_ac(ie,iw,ik) - vxcnn(ie,ik)
                sRe = dble(dsc)
                sIm = aimag(dsc)
                div = (w(iw)-enk-sRe)**2 + sIm**2
                sf(ie,iw,ik) = 1.d0/pi * abs(sIm) / div

            end do
            
        end do ! ie

    end do ! ik

    if (input%gw%selfenergy%actype == 'aaa') then
        call delete_aaa_approximant(aaa_plus)
        call delete_aaa_approximant(aaa_minus)
    end if

    ! Format specification
    n = nbgw-ibgw+1
    write(frmt1,'("(",i8,"f14.6)")') 1+n
    write(frmt2,'("(",i8,"f14.6)")') 1+2*n
    ! write(*,*) trim(frmt1), trim(frmt2)

    ! Output
    open(70,file='SpectralFunction.dat',form='FORMATTED',status='UNKNOWN',action='WRITE')
    open(71,file='SelfC-AC-Re.dat',form='FORMATTED',status='UNKNOWN',action='WRITE')
    open(72,file='SelfC-AC-Im.dat',form='FORMATTED',status='UNKNOWN',action='WRITE')
    open(73,file='Delta.dat',form='FORMATTED',status='UNKNOWN',action='WRITE')
    open(74,file='SelfC-Re.dat',form='FORMATTED',status='UNKNOWN',action='WRITE')
    open(75,file='SelfC-Im.dat',form='FORMATTED',status='UNKNOWN',action='WRITE')
    do ik = 1, kset%nkpt
        write(70,*) '# ik = ', ik
        write(71,*) '# ik = ', ik
        write(72,*) '# ik = ', ik
        write(73,*) '# ik = ', ik
        do iw = 1, nw
            write(70,trim(frmt1)) w(iw)+eferks, sf(:,iw,ik)
            write(71,trim(frmt1)) w(iw), dble(sc_ac(:,iw,ik))
            write(72,trim(frmt1)) w(iw), aimag(sc_ac(:,iw,ik))
            write(73,trim(frmt1)) w(iw), dble(w(iw)-evalks(:,ik)-selfex(:,ik)+vxcnn(:,ik))
        end do
        write(70,*); write(70,*)
        write(71,*); write(71,*)
        write(72,*); write(72,*)
        write(73,*); write(73,*)
        !---------------------------
        write(74,*) '# ik = ', ik
        write(75,*) '# ik = ', ik
        do iom = -freq_selfc%nomeg, freq_selfc%nomeg
            if (iom < 0) then
                write(74,trim(frmt1)) -freq_selfc%freqs(abs(iom)), dble(conjg(selfec(:,abs(iom),ik)))
                write(75,trim(frmt1)) -freq_selfc%freqs(abs(iom)), aimag(conjg(selfec(:,abs(iom),ik)))
            else if (iom > 0) then
                write(74,trim(frmt1)) freq_selfc%freqs(iom), dble(selfec(:,iom,ik))
                write(75,trim(frmt1)) freq_selfc%freqs(iom), aimag(selfec(:,iom,ik))
            end if
        end do
        write(74,*); write(74,*)
        write(75,*); write(75,*)
    end do
    close(70)
    close(71)
    close(72)
    close(73)
    close(74)

    ! clear memory
    deallocate(w, zj, fj)
    call delete_selfenergy
    deallocate(vxcnn)
    deallocate(sc_ac, sf)
      
    return
end subroutine