
subroutine calc_coulpot_mb

    use modmain
    use modgw
    use mod_hartreefock
    use modmpi

    implicit none
    integer(4) :: iq, i, j

!------------------------------------------------!
! Intialize Product Mixed Basis
!------------------------------------------------!
    call init_mixed_basis
    
!-----------------------------------------------------
! Calculate the integrals to treat the singularities at G+q->0
!-----------------------------------------------------
    call setsingc
!-----------------------------------------------------
! Precalculate Coulomb potential in MB representation
!-----------------------------------------------------
    if (rank == 0) then
        call boxmsg(fgw,'-','Precalculate Coulomb potential in MB representation')
    end if
    if (allocated(vcmbsiz)) deallocate(vcmbsiz)
    allocate(vcmbsiz(nqptnr))
    vcmbsiz(:) = 0
    if (allocated(vcmb)) deallocate(vcmb)
    allocate(vcmb(matsizmax,matsizmax,nqptnr))
    vcmb(:,:,:) = zzero
#ifdef MPI
    do iq = firstofset(rank,nqptnr), lastofset(rank,nqptnr)
#else
    do iq = 1, nqptnr
#endif
        Gamma = gammapoint(iq)
! Set the size of the basis for the corresponding q-point
        matsiz=locmatsiz+ngq(iq)
        if (rank==0) then
            write(fgw,101) iq, locmatsiz, ngq(iq), matsiz
        end if
! Calculate the interstitial mixed basis functions      
        call diagsgi(iq)
! Calculate the transformation matrix between pw's and the mixed basis functions
        call calcmpwipw(iq)
! Calculate the bare Coulomb matrix
        call calcbarcmb(iq)
! Reduce MB size by choosing eigenvectors of barc with eigenvalues larger than barcevtol
        call setbarcev(input%gw%BareCoul%barcevtol)
! store into global array
        vcmbsiz(iq) = mbsiz
        vcmb(1:matsiz,1:mbsiz,iq) = barcvm(1:matsiz,1:mbsiz)

        if ((debug).and.(rank==0)) then
            call linmsg(fgw,'-',' Vc_MB ')
            write(fgw,*)'for q-point ', iq
            write(fgw,*)'vcmbsiz=', vcmbsiz
            write(fgw,*)'eigenvalues='
            do i = 1, matsizmax, 20
                write(fgw,*) i, barcev(i)
            end do
            write(fgw,*)'SUM(vcmb)=', sum(vcmb(:,:,iq))
            call linmsg(fgw,'-','')
        end if
        
    end do ! iq
#ifdef MPI
    call mpi_allgatherv_ifc(nqptnr,1,ibuf=vcmbsiz)
    call mpi_allgatherv_ifc(nqptnr,matsizmax*matsizmax,zbuf=vcmb)
    call barrier
#endif
! Clean up memory (global arrays which come from modgw)
    if(allocated(barc))deallocate(barc)
    if(allocated(sqbarc))deallocate(sqbarc)
    if(allocated(barcev))deallocate(barcev)
    if(allocated(vmat))deallocate(vmat)
    if(allocated(vbas))deallocate(vbas)
    if(allocated(barcvm))deallocate(barcvm)
    if(allocated(sgi))deallocate(sgi)
    if(allocated(mpwipw))deallocate(mpwipw)
    
101 format(10x,/, &
    &       'Data for q-point:',i4,//,10x,'Mixed basis:',/,10x, &
    &       'Number of atomic basis functions:       ',i4,/,10x, &
    &       'Number of interstitial basis functions: ',i4,/,10x, &
    &       'Total number of basis functions:        ',i4,/)
    return
end subroutine
