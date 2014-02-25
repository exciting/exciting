
subroutine deltax
    use modinput
    use modmain
    use modmpi
    implicit none
! local variables
    integer :: is,ia,ias,ist,band0,band1
    integer :: ilm,irc,ik,ir
    integer :: ihomo, ilumo
    complex(8) :: zvx
! allocatable arrays
    complex(8), allocatable :: vnlvv(:)
    complex(8), allocatable :: delta(:,:)
    complex(8), allocatable :: apwalm(:,:,:,:)
    complex(8), allocatable :: evecfv(:,:)
    complex(8), allocatable :: evecsv(:,:)
    complex(8), allocatable :: wfmt(:,:,:,:,:)
    complex(8), allocatable :: wfir(:,:,:)
    complex(8), allocatable :: zrhomt(:,:,:)
    complex(8), allocatable :: zrhoir(:)
! external functions
    complex(8) zfinp,zfmtinp
    external zfinp,zfmtinp
!_______________________________________________________________________________
!
! test option: range of output states
    ihomo = 0
    ilumo = 1000
    do ik = 1, nkpt
        call getevalsv(vkl(:,ik),evalsv(:,ik))
        do ist = 1, nstsv
            if ((evalsv(ist,ik).le.efermi).and.(evalsv(ist+1,ik).gt.efermi)) then
                ihomo = max(ihomo,ist)
                ilumo = min(ilumo,ist+1)
                exit
            end if
        end do
    end do
!    band0 = max(1,ihomo-5)
!    band1 = min(nstsv,ihomo+5)

     band0 = 1
     band1 = nstsv 

!    if (rank==0) then
!        write(*,*) "efermi=", efermi
!        write(*,*) "ihomo=", ihomo
!        write(*,*) "ilumo=", ilumo
!    end if

!------------------
! k-point loop
!------------------

    allocate(delta(band0:band1,nkpt))
    call barrier

#ifdef MPI
    do ik = firstofset(rank,nkpt), lastofset(rank,nkpt)
#else    
    do ik = 1, nkpt
#endif

! get the eigenvalues/vectors from file for input k-point
        call getevalsv(vkl(:,ik),evalsv(:,ik))

        allocate(evecfv(nmatmax,nstfv))
        call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv)

        allocate(evecsv(nstsv,nstsv))
        call getevecsv(vkl(:,ik),evecsv)

! find the matching coefficients
        allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
        call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),apwalm)

! calculate the wavefunctions for all states
        allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
        allocate(wfir(ngrtot,nspinor,nstsv))
        call genwfsv(.false.,ngk(1,ik),igkig(:,1,ik),evalsv(:,ik), &
       &  apwalm,evecfv,evecsv,wfmt,wfir)
        deallocate(apwalm,evecfv,evecsv)

! 1 step:
! calculate <k|Vx_nl|k> matrix elements
        allocate(vnlvv(nstsv))        
        !call oepvnlk_diag(ik,vnlvv)
        call oepvnlk_full(ik,vnlvv)

! 2 step:
! Calculation of  <k|Vx|k>
        allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot))
        allocate(zrhoir(ngrtot))
        do ist = band0, band1
            do is = 1, nspecies
                do ia = 1, natoms(is)
                    ias = idxas(ia,is)
                    do ilm = 1, lmmaxvr
                        do irc = 1, nrcmtmax
                            zrhomt(ilm,irc,ias) = zvxmt(ilm,irc,ias)*wfmt(ilm,irc,ias,1,ist)
                        end do
                    end do
                 end do
             end do
             do ir = 1, ngrtot
                 zrhoir(ir) = zvxir(ir)*wfir(ir,1,ist)
             end do
             zvx = zfinp(.true.,wfmt(:,:,:,:,ist),zrhomt,wfir(:,:,ist),zrhoir)
             delta(ist,ik) = vnlvv(ist)-zvx  
! end loop over iband
        end do
        deallocate(vnlvv)
        deallocate(wfmt,wfir,zrhomt,zrhoir)

! end loop over ik
    end do

#ifdef MPI
    call mpi_allgatherv_ifc(nkpt, band1-band0+1, zbuf=delta)
    call barrier
#endif

    if (rank==0) then
          open(500,file='DELTAX'//trim(filext),action='WRITE',form='FORMATTED')
          write(500,*) nkpt, band1
          do ik=1,nkpt
              write(500,'(3f12.8,i5)') vkl(1,ik), vkl(2,ik), vkl(3,ik)
              do ist = band0, band1
                  write(500,'(6e20.6)') dble(delta(ist,ik))
              end do
              write(500,*)
          ! end loop over ik
          end do
          close(500)
          write(*,*)
          write(*,*) "Info(deltax): OEP exchange potential discontinuity is printed in DELTAX.OUT"
          write(*,*)
    end if

    deallocate(delta)

return
end subroutine

