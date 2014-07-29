
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
    complex(8) :: mvm, sum
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
        
    ! Calculate the integration weights using the linearized tetrahedron method
    if (allocated(evaldft)) deallocate(evaldft)
    allocate(evaldft(nstsv,nkpt))
    evaldft(:,:) = evalsv(:,:)
    call kintw
    deallocate(evaldft)
    
    ! Initialize mixed product basis
    call init_mb

! ***** revert back from GW definition of rwfcr
    do is=1,nspecies
        do ia=1,natoms(is)
            ias=idxas(ia,is)
            do ist=1,ncore(is)
                l=spl(ist,is)
                norm=sqrt(0.5d0*spocc(ist,is)/dble(2*l+1))
                do ir=1,nrmt(is)
                    rwfcr(ir,1,ist,ias)=spr(ir,is)*rwfcr(ir,1,ist,ias)/norm
                end do
            end do ! ist
        end do
    end do
    
!---------------------------------------
! Loop over k-points
!---------------------------------------

    exnl = 0.d0

#ifdef MPI
    do ikp = firstk(rank), lastk(rank)
#else
    do ikp = 1, nkpt
#endif    
        ik = idikp(ikp)
       
        vxnl(:,:,ikp) = zzero

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

!------------------------------------------------------------
! (Re)-Calculate the M^i_{nm}(k,q) matrix elements for given k and q
!------------------------------------------------------------
            call expand_basis(ik,iq)

!------------------------------------------------------------     
! (Re)-Calculate non-local (Hartree-Fock) potential
!------------------------------------------------------------

! Diagonalize the bare Coulomb matrix
            call setbarcev(0.d0)

! Valence electron contribution
            allocate(minm(mbsiz,nstfv,nstfv))
            nmdim = nstfv*nstfv
            call zgemm('c','n',mbsiz,nmdim,matsiz, &
            &          zone,barcvm,matsiz,minmmat,matsiz, &
            &          zzero,minm,mbsiz)

            jk = kqid(ik,iq)
            do ie1 = 1, nstsv
                do ie2 = 1, nstsv
                    sum = zzero
                    do ie3 = 1, nstsv
                        if (dabs(kiw(ie3,jk))>1.d-6) then
                            mvm = zdotc(mbsiz,minm(1:mbsiz,ie1,ie3),1,minm(1:mbsiz,ie2,ie3),1)
                            sum = sum+mvm*kiw(ie3,jk)
                        end if
                    end do ! ie3
                    vxnl(ie1,ie2,ikp) = vxnl(ie1,ie2,ikp)-sum
                end do ! ie2
                if (Gamma .and. (evalsv(ie1,ikp)<=efermi)) then
                    vxnl(ie1,ie1,ikp) = vxnl(ie1,ie1,ikp)-4.d0*pi*vi*singc2
                end if
            end do ! ie1
                
            deallocate(minm)
            deallocate(minmmat)

! Core electron contribution
            if (iopcore<=1) then
                
                allocate(minm(mbsiz,nstfv,ncg))
                nmdim = nstfv*ncg
                call zgemm('c','n',mbsiz,nmdim,locmatsiz, &
                &          zone,barcvm,matsiz,mincmat,locmatsiz, &
                &          zzero,minm,mbsiz)
                
                do ie1 = 1, nstsv
                    do ie2 = 1, nstsv
                        sum = zzero
                        do icg = 1, ncg
                            mvm = zdotc(mbsiz,minm(1:mbsiz,ie1,icg),1,minm(1:mbsiz,ie2,icg),1)
                            ! BZ integration weight
                            is = corind(icg,1)
                            ia = corind(icg,2)
                            ias= idxas(ia,is)
                            ic = corind(icg,3)
                            sum = sum+mvm*ciw(ias,ic)
                        enddo ! ie3
                        vxnl(ie1,ie2,ikp) = vxnl(ie1,ie2,ikp)-sum
                    end do ! ie2
                end do ! ie1
                    
                deallocate(minm)
                deallocate(mincmat)
                    
            end if ! iopcore

        end do ! iq

!------------------------------------------------------------
! Evaluate HF energy E_HF(k)
!------------------------------------------------------------
        do ie1 = 1, nstsv
            if (evalsv(ie1,ikp)<=efermi) then
                exnl = exnl+occmax*wkpt(ikp)*vxnl(ie1,ie1,ikp)*nkptnr
            end if
        end do

!------------------------------------------------------------
! Debugging Info
!------------------------------------------------------------
        if ((debug).and.(rank==0)) then
            call linmsg(fgw,'-','')
            call linmsg(fgw,'-',' Diagonal elements of Vx_NL_nn ')
            write(fgw,*) 'for k-point ', ikp
            do ie1 = 1, nstfv
                write(fgw,*) ie1, vxnl(ie1,ie1,ikp)
            end do
            call linmsg(fgw,'-','')
            call linmsg(fgw,'-',' Vx_NL_nn ')
            do ie1 = 1, nstfv
                write(fgw,*) vxnl(ie1,:,ikp)
            end do
            call linmsg(fgw,'-','')
        end if
    
    end do ! ikp

#ifdef MPI
    call MPI_ALLREDUCE(MPI_IN_PLACE, exnl, 1, MPI_DOUBLE_PRECISION, &
    &                  MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, vxnl, nstfv*nstfv*nkpt, MPI_DOUBLE_COMPLEX, &
    &                  MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

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
