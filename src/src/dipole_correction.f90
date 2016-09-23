subroutine dipole_correction(vdplmt,vdplir)

    use modinput
    use modmain
    use modmpi
    
    implicit none
    real(8), intent(out) :: vdplmt(lmmaxvr,nrmtmax,natmtot)
    real(8), intent(out) :: vdplir(ngrtot)
    
    integer :: is, ia, ias, ir
    integer :: i1, i2, i3
    real(8) :: dplmt, dplir, dpl, zm, z, z0
    real(8) :: t1, t2, rv(3), norm(3), A
    real(8) :: fr(nrmtmax), gr(nrmtmax), cf(3,nrmtmax)
    ! external function
    real(8) :: rfinp
    integer :: ivac
    
    integer, allocatable :: idx(:)
    real(8), allocatable :: tmat(:,:,:), rhoz(:)

    ! vacuum region direction x/y/z -> 1/2/3
    ivac = 3

    !--------------------------
    ! Find the surface area A
    !--------------------------
    select case (ivac)
      case(1)
        call r3cross(input%structure%crystal%basevect(:,2), &
        &            input%structure%crystal%basevect(:,3), &
        &            norm(:))
      case(2)
        call r3cross(input%structure%crystal%basevect(:,3), &
        &            input%structure%crystal%basevect(:,1), &
        &            norm(:))
      case(3)
        call r3cross(input%structure%crystal%basevect(:,1), &
        &            input%structure%crystal%basevect(:,2), &
        &            norm(:))
    end select

    A = dsqrt(dot_product(norm(:),norm(:)))
    if (rank==0) write(*,*) 'Surface area A=', A
    
    zm = dsqrt(dot_product(input%structure%crystal%basevect(:,ivac), &
    &                      input%structure%crystal%basevect(:,ivac)))
    if (rank==0) write(*,*) 'zm=', zm
    
    ! position of the dipole potential discontinuity
    z0 = input%groundstate%dipoleposition*zm
    
    t1 = 1000.d0
    t2 = 0.d0
    do i3 = 0, ngrid(3)-1
    do i2 = 0, ngrid(2)-1
    do i1 = 0, ngrid(1)-1
      rv(:) = dble(i1)/dble(ngrid(1))*input%structure%crystal%basevect(:,1)+ &
      &       dble(i2)/dble(ngrid(2))*input%structure%crystal%basevect(:,2)+ &
      &       dble(i3)/dble(ngrid(3))*input%structure%crystal%basevect(:,3)
      z = rv(ivac)
      if (abs(z-z0)<t1) then
        t1 = abs(z-z0)
        t2 = z
      end if
    end do
    end do
    end do
    if (rank==0) then
      write(*,*) 'Starting position: z0=', z0
      write(*,*) 'Nearest FFT grid point', t2
    end if
    z0 = t2 !!!!
    
if (.false.) then 
    ! determine minimum of the plane averaged IR density
    allocate(tmat(ngrid(1),ngrid(2),ngrid(3)))
    ir = 0
    do i3 = 1, ngrid(3)
    do i2 = 1, ngrid(2)
    do i1 = 1, ngrid(1)
      ir = ir+1
      tmat(i1,i2,i3) = rhoir(ir) 
    end do
    end do
    end do
    ! plane averaging
    allocate(rhoz(ngrid(3)))
    t1 = 1.d0/dble(ngrid(1)*ngrid(2))
    do i3 = 1, ngrid(3)
      rhoz(i3) = sum(tmat(:,:,i3))*t1
      if (rank==0) write(77,*) i3, rhoz(i3)
    end do
    deallocate(tmat,rhoz)
end if
    
    !==============================
    ! Calculate the dipole moment
    !==============================

    ! MT part
    dplmt = 0d0
    do is = 1, nspecies
    do ia = 1, natoms(is)
      ias = idxas(ia,is)
      z = atposc(ivac,ia,is)
      if (z>z0) z = z-zm
      ! 00
      do ir = 1, nrmt(is)
        fr(ir) = rhomt(1,ir,ias)*spr(ir,is)**2
      end do ! ir
      call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
      t1 = (dsqrt(fourpi)*gr(nrmt(is))+spzn(is))*z
      dplmt = dplmt+t1
      ! 10
      do ir = 1, nrmt(is)
        fr(ir) = rhomt(3,ir,ias)*spr(ir,is)**3
      end do ! ir
      call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
      t1 = dsqrt(fourpi/3d0)*gr(nrmt(is))
      dplmt = dplmt+t1
    end do ! ia
    end do ! is
    if (rank==0)  write(*,*) 'dplmt = ', dplmt/A
   
    ! IR part
    dplir = 0.d0
    ir = 0
    do i3 = 0, ngrid(3)-1
    do i2 = 0, ngrid(2)-1
    do i1 = 0, ngrid(1)-1
      ir = ir+1
      rv(:) = dble(i1)/dble(ngrid(1))*input%structure%crystal%basevect(:,1)+ &
      &       dble(i2)/dble(ngrid(2))*input%structure%crystal%basevect(:,2)+ &
      &       dble(i3)/dble(ngrid(3))*input%structure%crystal%basevect(:,3)
      z = rv(ivac)
      if (abs(z-z0)<1.d-8) then
        z = 0
      else if (z>z0) then
        z = z-zm
      end if
      dplir = dplir+rhoir(ir)*dble(cfunir(ir))*z
    end do
    end do
    end do
    dplir = dplir*omega/dble(ngrtot)
    if (rank==0) write(*,*) 'dplir = ', dplir/A
    
    dpl = (dplmt+dplir)/A
    if (rank==0) then
      write(60,*)
      write(60,'(" Slab surface-dipole density", T45, ": ", F18.8)') dpl
      write(*,'(" Slab surface-dipole density", T45, ": ", F18.8)') dpl
      write(60,*)
    end if

    !===============================================
    ! Dipole correction for the Coulomb potential
    !===============================================
    
    vdplmt(:,:,:) = 0.d0
    vdplir(:) = 0.d0
    
    ! constant prefactor
    t1 = fourpi*dpl
    
    ! MT part
    do is = 1, nspecies
    do ia = 1, natoms(is)
      ias = idxas(ia,is)
      z = atposc(ivac,ia,is)
      if (z>z0) z = z-zm
      do ir = 1, nrmt(is)
        ! 00
        vdplmt(1,ir,ias) = t1*dsqrt(fourpi)*(z/zm-0.5d0)
        ! 10
        vdplmt(3,ir,ias) = t1*dsqrt(fourpi/3d0)*spr(ir,is)/zm
      end do ! ir
    end do ! ia
    end do ! is
    
    ! IR part
    ir = 0
    do i3 = 0, ngrid(3)-1
    do i2 = 0, ngrid(2)-1
    do i1 = 0, ngrid(1)-1
      ir = ir+1
      rv(:) = dble(i1)/dble(ngrid(1))*input%structure%crystal%basevect(:,1)+ &
      &       dble(i2)/dble(ngrid(2))*input%structure%crystal%basevect(:,2)+ &
      &       dble(i3)/dble(ngrid(3))*input%structure%crystal%basevect(:,3)
      z = rv(ivac)
      if (abs(z-z0)<1.d-8) then
        z = 0
      else if (z>z0) then
        z = z-zm
      end if
      vdplir(ir) =  t1*(z/zm-0.5d0)
    end do
    end do
    end do
    
    !================================================
    ! Calculate the dipole correction energy Eq.(13)
    !================================================
    
    ! ionic contribution
    endipc = 0.d0
    do is = 1, nspecies
    do ia = 1, natoms(is)
      ias = idxas(ia,is)
      z = atposc(ivac,ia,is)
      endipc = endipc+0.5d0*spzn(is)*fourpi*dpl*(z/zm-0.5d0)
    end do ! ia
    end do ! is
    if (rank==0) write(*,*) 'Ionic: = ', endipc
    
    ! electronic contribution
    t1 = 0.5d0*rfinp(1,rhomt,vdplmt,rhoir,vdplir)
    if (rank==0) write(*,*) 'Electronic: = ', t1
    
    endipc = endipc-t1
    if (rank==0) write(*,*) 'Dipole correction energy = ', endipc    
    
    return
end subroutine

