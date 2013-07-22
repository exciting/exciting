
subroutine calc_vxnl
    
    use modmain
    use modgw
    use modfvsystem
    use mod_hartreefock
    use modmpi
            
    implicit none

    integer(4) :: ikp, ik, iq
    integer(4) :: ie1, ie2, ie3, icg
    real(8)    :: tstart, tend
    complex(8) :: t1, mvm
    
    real(8),    allocatable :: eval(:)
    complex(8), allocatable :: vnl(:,:)
    
    complex(8), external :: zdotc

!------------------------------------------------!
! LAPW basis initialization
!------------------------------------------------!
! generate the core wavefunctions and densities
    Call gencore
! compute linearisation energies
    Call linengy
    if (rank.eq.0) Call writelinen
! generate the APW radial functions
    Call genapwfr
! generate the local-orbital radial functions
    Call genlofr(.false.)
! compute the overlap radial integrals
    Call olprad
! compute the Hamiltonian radial integrals
    Call hmlrad

!------------------------------------------------!
!   Calculate/Update Coulomb Potential in MB representation
!------------------------------------------------!
    if (rank .Eq. 0) call boxmsg(60,'-','Coulomb potential in MB')
    call cpu_time(tstart)
    call calc_coulpot_mb
    call cpu_time(tend)
    if (rank == 0) then
        call write_cputime(60,tend-tstart,'calc_coul_mb')
        call boxmsg(60,'-','')
    end if

!------------------------------------------------!
!   Calculate/Update the non-local exchange potential
!------------------------------------------------!
    if (rank == 0) then
        call boxmsg(fgw,'-','Calculate Vx_NL')
    end if
    
    allocate(eval(nstsv))
    call readfermi
    nomax = 0
    do ikp = 1, nkpt
        call getevalsv(vkl(:,ikp),eval)
        do ie1 = 1, nstsv
            if (eval(ie1)<=efermi) nomax=max(ie1,nomax)
        end do
    end do

!---------------------------------------
! Loop over k-points
!---------------------------------------
    do ikp = 1, nkpt
        ik = idikp(ikp)

        if (allocated(vnl)) deallocate(vnl)
        allocate(vnl(nstfv,nstfv))
        vnl(:,:)=zzero

!---------------------------------------
! Integration over BZ
!---------------------------------------
        do iq = 1, nqptnr

            Gamma = gammapoint(iq)
! Set the size of the basis for the corresponding q-point
            matsiz=locmatsiz+ngq(iq)
            if (rank == 0) write(fgw,101) iq, locmatsiz, ngq(iq), matsiz

! Calculate the interstitial mixed basis functions
            call diagsgi(iq)
! Calculate the transformation matrix between pw's and the mixed basis functions
            call calcmpwipw(iq)

!------------------------------------------------------------
! calculate the M^i_{nm}(k,q) matrix elements for given k and q
!------------------------------------------------------------
            call expand_basis(ik,iq)

!------------------------------------------------------------     
! Calculate non-local (Hartree-Fock) potential
!------------------------------------------------------------     
! Valence electron contribution
            do ie1 = 1, nstsv
                do ie2 = 1, nstsv
                    t1 = zzero
                    do ie3 = 1, nomax
                        mvm = zdotc(mbsiz,minmmat(1:mbsiz,ie1,ie3),1,minmmat(1:mbsiz,ie2,ie3),1)
                        t1 = t1+mvm
                    enddo ! ie3
                    vnl(ie1,ie2) = vnl(ie1,ie2)-t1/dble(nqptnr)
                end do ! ie2
                if (Gamma.and.(ie1<=nomax)) then
                    vnl(ie1,ie1) = vnl(ie1,ie1)-4.d0*pi*vi*singc2
                end if
            end do ! ie1
            deallocate(minmmat)

! Core electron contribution
            if (iopcore <= 1) then
                do ie1 = 1, nstsv
                    do ie2 = 1, nstsv
                        t1 = zzero
                        do icg = 1, ncg
                            mvm = zdotc(mbsiz,mincmat(1:mbsiz,ie1,icg),1,mincmat(1:mbsiz,ie2,icg),1)
                            t1 = t1+mvm
                        enddo ! ie3
                        vnl(ie1,ie2) = vnl(ie1,ie2)-t1/dble(nqptnr)
                    end do ! ie2
                end do ! ie1
                deallocate(mincmat)
            end if ! iopcore
          
        end do ! iq

!------------------------------------------------------------
! Store k-dependent non-local potential into global array
!------------------------------------------------------------
        vxnl(:,:,ikp) = vnl(:,:)
    
!------------------------------------------------------------
! Evaluate HF energy E_HF(k)
!------------------------------------------------------------
        exnlk(ikp) = 0.d0
        do ie1 = 1, nomax
            exnlk(ikp) = exnlk(ikp)+0.5d0*vnl(ie1,ie1)*occmax
        end do

!------------------------------------------------------------
! Debugging Info
!------------------------------------------------------------
        if ((debug).and.(rank==0)) then
            write(fgw,*)
            write(fgw,*) 'Hartree-Fock energy for ikp=', ikp
            write(fgw,*) 'exnl=', exnlk(ikp)
            call linmsg(fgw,'-','')
            call linmsg(fgw,'-',' Diagonal elements of Vx_NL_nn ')
            write(fgw,*) 'for k-point ', ikp
            do ie1 = 1, nstfv
                write(fgw,*) ie1, vnl(ie1,ie1)
            end do
            call linmsg(fgw,'-','')
            call linmsg(fgw,'-',' Vx_NL_nn ')
            do ie1 = 1, nstfv
                write(fgw,*) vnl(ie1,:)
            end do
            call linmsg(fgw,'-','')
        end if
    
    end do ! ikp

101 format(10x,/, &
    &       'Data for q-point:',i4,//,10x,'Mixed basis:',/,10x, &
    &       'Number of atomic basis functions:       ',i4,/,10x, &
    &       'Number of interstitial basis functions: ',i4,/,10x, &
    &       'Total number of basis functions:        ',i4,/)
      
    return
end subroutine
