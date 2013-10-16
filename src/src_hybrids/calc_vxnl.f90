
subroutine calc_vxnl
    
    use modmain
    use modgw
    use modfvsystem
    use mod_hybrids
    use modmpi
            
    implicit none

    integer(4) :: ikp, ik, iq, ikq
    integer(4) :: ie1, ie2, ie3, icg
    integer(4) :: ia, is, ias, ir, ist, l
    real(8)    :: norm
    real(8)    :: tstart, tend
    complex(8) :: t1, mvm
    real(8), allocatable :: eval(:)    
    complex(8), allocatable :: vnl(:,:)
    
    complex(8), external :: zdotc

    call cpu_time(tstart)

!------------------------------------------------!
!   Calculate/Update the non-local exchange potential
!------------------------------------------------!
    if (rank == 0) then
        call boxmsg(fgw,'-','Calculate Vx_NL')
    end if

! State-index of VBM    
    allocate(eval(nstsv))
    call readfermi
    nomax = 0
    do ikp = 1, nkpt
        call getevalsv(vkl(:,ikp),eval)
        do ie1 = 1, nstsv
            if (eval(ie1)>efermi) then
                nomax = max(ie1-1,nomax)
                exit ! ie1 loop
            end if
        end do
    end do
    deallocate(eval)
    
!---------------------------------------
!   Mixed basis initialization (GW routine)
!---------------------------------------
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
    call barrier

    ikq = 0
    do ikp = 1, nkpt
        ik = idikp(ikp)

        if (allocated(vnl)) deallocate(vnl)
        allocate(vnl(nstfv,nstfv))
        vnl(:,:)=zzero

!---------------------------------------
! Integration over BZ
!---------------------------------------
        do iq = 1, nqptnr

! decide if point is done by this process
            ikq = ikq+1
            if (mod(ikq-1,procs) == rank) then

                
! Set the size of the basis for the corresponding q-point
                matsiz=locmatsiz+ngq(iq)
                if (rank==0) write(fgw,101) ikp, iq, locmatsiz, ngq(iq), matsiz

!------------------------------------------------------------
! (Re)-Calculate the Coulomb matrix elements in MB representation
!------------------------------------------------------------
                Gamma = gammapoint(iq)

! Calculate the interstitial mixed basis functions
                call diagsgi(iq)

! Calculate the transformation matrix between pw's and the mixed basis functions
                call calcmpwipw(iq)

! Calculate the bare Coulomb matrix
                call calcbarcmb(iq)

! Reduce MB size by choosing eigenvectors of barc with eigenvalues larger than barcevtol
                call setbarcev(input%gw%BareCoul%barcevtol)            

!------------------------------------------------------------
! (Re)-Calculate the M^i_{nm}(k,q) matrix elements for given k and q
!------------------------------------------------------------
                call expand_basis(ik,iq)

!------------------------------------------------------------     
! (Re)-Calculate non-local (Hartree-Fock) potential
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
                    if (Gamma .and. (ie1 <= nomax)) then
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

            end if !   
        end do ! iq

!------------------------------------------------------------
! Store k-dependent non-local potential into global array
!------------------------------------------------------------
        vxnl(:,:,ikp) = vnl(:,:)
        deallocate(vnl)
    
!------------------------------------------------------------
! Evaluate HF energy E_HF(k)
!------------------------------------------------------------
        exnlk(ikp) = 0.d0
        do ie1 = 1, nomax
            exnlk(ikp) = exnlk(ikp)+0.5d0*vxnl(ie1,ie1,ikp)*occmax
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
    call MPI_ALLREDUCE(MPI_IN_PLACE, exnlk, nkpt, MPI_DOUBLE_PRECISION, &
   &  MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, vxnl, nstfv*nstfv*nkpt, MPI_DOUBLE_COMPLEX, &
   &  MPI_SUM, MPI_COMM_WORLD, ierr)
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
