
subroutine calc_vxnl
    
    use modmain
    use modgw
    use modfvsystem
    use mod_hybrids
    use modmpi
            
    implicit none

    integer(4) :: ikp, ik, jk, iq
    integer(4) :: ie1, ie2, ie3, icg, ic
    integer(4) :: ia, is, ias, ir, ist, l
    integer(4) :: n, nmdim
    real(8)    :: norm, egap
    real(8)    :: tstart, tend
    complex(8) :: mvm, sum, fnk
    complex(8), allocatable :: minm(:,:,:)
    complex(8), allocatable :: vnl(:,:)
    
    complex(8), external :: zdotc

    call cpu_time(tstart)
    if (rank == 0) call boxmsg(fgw,'-','Calculate Vx_NL')
    
    !------------------------------------------!
    ! Matrix elements of non-local potential   !
    !------------------------------------------!
    if (allocated(vxnl)) deallocate(vxnl)
    allocate(vxnl(nstfv,nstfv,nkpt))
    vxnl(:,:,:) = zzero
        
    ! Calculate the integration weights using the linearized tetrahedron method
    !if (allocated(evaldft)) deallocate(evaldft)
    !allocate(evaldft(nstsv,nkpt))
    !evaldft(:,:) = evalsv(:,:)
    !call kintw
    !deallocate(evaldft)
    !deallocate(kcw)
    !if (iopcore<2) deallocate(unw)
    
    ! Initialize mixed product basis
    call init_mb
    
!---------------------------------------
! Loop over k-points
!---------------------------------------

#ifdef MPI
    do ikp = firstk(rank), lastk(rank)
#else
    do ikp = 1, nkpt
#endif    
        ik = idikp(ikp)
       
!---------------------------------------
! Integration over BZ
!---------------------------------------
        do iq = 1, nqptnr
          
            Gamma = gammapoint(iq)

! Set the size of the basis for the corresponding q-point
            matsiz = locmatsiz+ngq(iq)
            if (rank==0) write(fgw,101) ikp, iq, locmatsiz, ngq(iq), matsiz

! Calculate the interstitial mixed basis functions
            call diagsgi(iq)

! Calculate the transformation matrix between pw's and the pw mixed basis functions
            call calcmpwipw(iq)

!------------------------------------               
! Calculate the bare Coulomb matrix
!------------------------------------
            call calcbarcmb(iq)
            call setbarcev(0.d0)

!------------------------------------------------------------
! (Re)-Calculate the M^i_{nm}(k,q) matrix elements for given k and q
!------------------------------------------------------------
            call expand_basis(ik,iq)

!------------------------------------------------------------     
! (Re)-Calculate non-local (Hartree-Fock) potential
!------------------------------------------------------------

! Valence electron contribution
            allocate(minm(mbsiz,nstfv,nstfv))
            nmdim = nstfv*nstfv
            call zgemm('c','n',mbsiz,nmdim,matsiz, &
            &          zone,barcvm,matsiz,minmmat,matsiz, &
            &          zzero,minm,mbsiz)
            deallocate(minmmat)

            jk = kqid(ik,iq)
            do ie1 = 1, nstsv
              do ie2 = ie1, nstsv
                sum = zzero
                do ie3 = 1, nstsv
                  if (evalsv(ie3,ikp)<=efermi) then
                    mvm = zdotc(mbsiz,minm(1:mbsiz,ie1,ie3),1,minm(1:mbsiz,ie2,ie3),1)
                    sum = sum+mvm
                  end if
                end do ! ie3
                vxnl(ie1,ie2,ikp) = vxnl(ie1,ie2,ikp)-sum/dble(nqptnr)
              end do ! ie2
            end do ! ie1
            deallocate(minm)

! Core electron contribution
            if (iopcore<=1) then
                
              allocate(minm(mbsiz,nstfv,ncg))
              nmdim = nstfv*ncg
              call zgemm('c','n',mbsiz,nmdim,locmatsiz, &
              &          zone,barcvm,matsiz,mincmat,locmatsiz, &
              &          zzero,minm,mbsiz)
              deallocate(mincmat)
                
              do ie1 = 1, nstsv
                do ie2 = ie1, nstsv
                  sum = zzero
                  do icg = 1, ncg
                    mvm = zdotc(mbsiz,minm(1:mbsiz,ie1,icg),1,minm(1:mbsiz,ie2,icg),1)
                    sum = sum+mvm
                  enddo ! ie3
                  vxnl(ie1,ie2,ikp) = vxnl(ie1,ie2,ikp)-sum/dble(nqptnr)
                end do ! ie2
              end do ! ie1
              deallocate(minm)
                    
            end if ! iopcore

        end do ! iq
        
        do ie1 = 1, nstsv
          do ie2 = ie1+1, nstsv
            vxnl(ie2,ie1,ikp) = conjg(vxnl(ie1,ie2,ikp))
          end do
        end do

!------------------------------------------------------------
! Debugging Info
!------------------------------------------------------------
        if ((.true.).and.(rank==0)) then
            call linmsg(fgw,'-','')
            call linmsg(fgw,'-',' Diagonal elements of Vx_NL_nn ')
            write(fgw,*) 'for k-point ', ikp
            do ie1 = 1, nstfv
                write(fgw,*) ie1, vxnl(ie1,ie1,ikp)
            end do
            call linmsg(fgw,'-','')
            !call linmsg(fgw,'-',' Vx_NL_nn ')
            !do ie1 = 1, nstfv
            !    write(fgw,*) vxnl(ie1,:,ikp)
            !end do
            !call linmsg(fgw,'-','')
        end if
    
    end do ! ikp

#ifdef MPI
    call mpi_allgatherv_ifc(nkpt,nstfv*nstfv,zbuf=vxnl)
    call barrier
#endif


    !-----------------------------------------
    ! Calculate the non-local exchange energy
    !-----------------------------------------
    exnl = 0.d0
    do ikp = 1, nkpt
      do ie1 = 1, nstfv
        if (evalsv(ie1,ikp)<=efermi) then
          fnk = occsv(ie1,ikp)/occmax
          vxnl(ie1,ie1,ikp) = vxnl(ie1,ie1,ikp)- &
          &                   4.d0*pi*vi*singc2* &
          &                   fnk
          exnl = exnl+occmax*fnk*wkpt(ikp)*vxnl(ie1,ie1,ikp)
        end if
      end do
    end do

101 format(10x,/, &
    &       'Data for (k,q)-point:',2i4,//,10x,                  &
    &       'Mixed basis:',/,10x,                                &
    &       'Number of atomic basis functions:       ',i4,/,10x, &
    &       'Number of interstitial basis functions: ',i4,/,10x, &
    &       'Total number of basis functions:        ',i4,/)

    call cpu_time(tend)
    if (rank==0) call write_cputime(fgw,tend-tstart, 'CALC_VXNL')
    
    return
end subroutine
