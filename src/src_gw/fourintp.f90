!BOP
!
!!ROUTINE: fourintp
!
!!INTERFACE:
!
subroutine fourintp(f1,nk1,kvecs1,f2,nk2,kvecs2,nb)
!      
! !DESCRIPTION:
!
! This subroutine interpolate function $f1_n( k )$ defined on the kmesh 1, to kmesh 2
! using 3D Smooth Fourier transform according to 
! PRB 38, 2721 (1988).
!
!!USES:
    use modinput
    use modmain, only: nsymcrys, pi, zzero, zone
    use fouri

!!LOCAL VARIABLES:
    implicit none
    integer(4), intent(in) :: nk1,nk2,nb
    real(8),    intent(in) :: kvecs1(3,nk1),kvecs2(3,nk2)
    complex(8), intent(in) :: f1(nk1,1:nb)
    complex(8), intent(out):: f2(nk2,1:nb) 
      
    integer(4) :: i, ist, ib, ik, jk, ir
    integer(4) :: info
    integer(4), allocatable  :: ipiv(:)
      
    real(8) :: den, pref, kdotr
    real(8) :: rmin,rlen,x2,x6,c1,c2
    real(8), dimension(3) :: r,rvec,kvec
    real(8), allocatable  :: rho(:)
      
    complex(8) :: expkr
    complex(8), allocatable :: dele(:,:)
    complex(8), allocatable :: h(:,:)
    complex(8), allocatable :: coef(:,:)
    complex(8), allocatable :: smat1(:,:),smat2(:,:)
    complex(8), allocatable :: sm2(:,:)

!!REVISION HISTORY:
!
! Created 02.03.06 by RGA
! Revisited July 2011 by DIN
!
!EOP
!BOC

    ! shortcut for basis vectors 
    rbas(:,1) = input%structure%crystal%basevect(:,1)
    rbas(:,2) = input%structure%crystal%basevect(:,2)
    rbas(:,3) = input%structure%crystal%basevect(:,3)

    ! Set rindex
    if(.not.setrindex_done) then  
      nsymcrys = 1
      call setrindex
      setrindex_done = .true.
    endif

    c1 = 0.25d0
    c2 = 0.25d0
    allocate(smat1(nk1,nst), &
    &        smat2(nk2,nst), &
    &        rho(nst),       &
    &        coef(nst,nb),   &
    &        ipiv(1:nk1-1),  &
    &        sm2(1:nk1-1,1:nst), &
    &        h(1:nk1-1,1:nk1-1), &
    &        dele(1:nk1-1,1:nb))

    den = dble(nsymcrys)

    ! Calculate the star expansion function at each irreducible k-point
    smat1(1:nk1,1:nst) = zzero
    do ik = 1, nk1
      kvec(1:3) = kvecs1(1:3,ik)
      do ir = 2, nrr
        ist = rst(1,ir)
        pref = dble(rst(2,ir))
        r(1:3) = dble(rindex(1:3,ir))
        kdotr = 2.0d0*pi*sum(r(1:3)*kvec(1:3)) 
        expkr = cmplx(cos(kdotr),sin(kdotr),8)
        smat1(ik,ist) = smat1(ik,ist)+pref*expkr/den
      enddo
    enddo 

    ! Calculate the curvature function (rho) for each star
    rho(1:nst) = 0.0d0
    ist = 1
    do ir = 2, nrr
      if (rst(1,ir).ne.ist) then
        ist = rst(1,ir)
        r(1:3) = dble(rindex(1:3,ir))
    	do i = 1, 3
    	  rvec(i) = r(1)*rbas(i,1)+r(2)*rbas(i,2)+r(3)*rbas(i,3)
        enddo
        rlen = sum(rvec(1:3)*rvec(1:3))
        if (ist.eq.2) rmin = rlen
        x2 = rlen/rmin
        x6 = x2*x2*x2
        rho(ist) = (1.0d0-c1*x2)*(1.0d0-c1*x2)+c2*x6
      endif
    enddo

    ! Set sm2(k)=smat(k)-smat(k_nkp) and dele
    do ik = 1, nk1-1
      do ist = 2, nst
        sm2(ik,ist) = smat1(ik,ist)-smat1(nk1,ist)
      enddo
      do ib = 1, nb
        dele(ik,ib) = f1(ik,ib)-f1(nk1,ib)
      enddo
    enddo

    ! Calculate the matrix H      
    h(1:nk1-1,1:nk1-1) = zzero
    do ik = 1, nk1-1
      do jk = 1, nk1-1
        do ist = 2, nst
          h(ik,jk) = h(ik,jk)+sm2(ik,ist)*conjg(sm2(jk,ist))/rho(ist)
        enddo
      enddo
    enddo

    ! Solve the Linear equations for the Lagrange multipliers
    call zgetrf(nk1-1,nk1-1,h,nk1-1,ipiv,info)
    call errmsg(info.ne.0,"fourintp","error when calling zgetrf")

    call zgetrs('n',nk1-1,nb,h,nk1-1,ipiv,dele,nk1-1,info)
    call errmsg(info.ne.0,"fourintp","error when calling zgetrs")

    ! Calculate the coefficients of the Star expansion
    coef(1,1:nb) = f1(nk1,1:nb)
    do ist = 2, nst
      coef(ist,1:nb) = zzero
      do ik = 1, nk1-1
        coef(ist,1:nb) = coef(ist,1:nb)+dele(ik,1:nb)*conjg(sm2(ik,ist))
      enddo
      coef(ist,1:nb) = coef(ist,1:nb)/rho(ist)
      coef(1,1:nb) = coef(1,1:nb)-coef(ist,1:nb)*smat1(nk1,ist)
    enddo

    ! Perform an interpolation to the new k-mesh
    smat2(1:nk2,1:nst) = zzero
    do ik = 1, nk2
      kvec(1:3) = kvecs2(1:3,ik)
      do ir = 1, nrr
        ist = rst(1,ir)
        pref = dble(rst(2,ir))
        r(1:3) = dble(rindex(1:3,ir))
        kdotr = 2.0d0*pi*sum(r(1:3)*kvec(1:3)) 
        expkr = cmplx(cos(kdotr),-sin(kdotr),8)
        smat2(ik,ist) = smat2(ik,ist)+pref*expkr/den
      enddo 
    enddo 
    call zgemm('n','n',nk2,nb,nst,zone,smat2,nk2,coef,nst,zzero,f2,nk2)
      
    deallocate(smat1,smat2,coef,rho,ipiv,sm2,h,dele)

    return
end subroutine
!EOC   
      
      
