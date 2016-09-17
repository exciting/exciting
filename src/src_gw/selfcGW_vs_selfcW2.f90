
subroutine selfcGW_vs_selfcW2

    use modinput
    use modmain, only : efermi, zzero
    use modgw
    use m_getunit
    implicit none
    
    integer :: ie, ikp, fid
    real(8) :: enk
    complex(8) :: ein, dsig, sigma
    complex(8) :: selfcGW, selfcW2
    ! parameters of the functions fitting the selfenergy
    integer :: npar
    complex(8), allocatable :: a(:), sc(:), poles(:)
    complex(8), allocatable :: sacpar(:,:,:)

    call getunit(fid)
    open(fid,file='SELFC-GW-W2.DAT',action='WRITE',form='FORMATTED')
    write(fid,*) '# Analytical continuation'
    write(fid,*) '# ik    ie    selfcGW    selfcW2    selfcGW-selfcW2'

    ! Parameters of perform the analytic function of the correlation self-energy
    npar = 2*input%gw%selfenergy%npol
    allocate(a(npar),sc(freq%nomeg),poles(npar))
    allocate(sacpar(npar,ibgw:nbgw,kset%nkpt))
    
    do ikp = 1, kset%nkpt
      do ie = ibgw, nbgw
      
        enk  = evalks(ie,ikp)-efermi
        ein = cmplx(enk,0.0d0,8)
        
        !---------------------------------------------
        ! GW self-energy
        !---------------------------------------------
        sc(1:freq%nomeg) = selfec(ie,1:freq%nomeg,ikp)
        call setsac(iopac,freq%nomeg, &
        &           npar,enk,sc,freq%freqs,a,poles)
        sacpar(:,ie,ikp) = a
        call getsac(iopac,freq%nomeg, &
        &           npar,enk,ein,freq%freqs,a,sigma,dsig)
        selfcGW = sigma
        
        !---------------------------------------------
        ! GW self-energy
        !---------------------------------------------
        sc(1:freq%nomeg) = selfecw2(ie,1:freq%nomeg,ikp)
        call setsac(iopac,freq%nomeg, &
        &           npar,enk,sc,freq%freqs,a,poles)
        sacpar(:,ie,ikp) = a
        call getsac(iopac,freq%nomeg, &
        &           npar,enk,ein,freq%freqs,a,sigma,dsig)
        selfcW2 = sigma
    
        write(fid,'(2i4,6f14.6)') ikp, ie, selfcGW, selfcW2, selfcGW-selfcW2
      
      enddo ! ie
      write(fid,*)
    enddo ! ikp
    
    if (allocated(a)) deallocate(a)
    if (allocated(sc)) deallocate(sc)
    if (allocated(poles)) deallocate(poles)
    if (allocated(sacpar)) deallocate(sacpar)
    
    close(fid)
    
    return
end subroutine
          
