
module mod_rgrid
  
  use modinput
  implicit none

  type rgrid
    integer :: npt
    integer :: iflag                     ! grid type flag: <0 - general periodic (xcrysden format)
    integer, allocatable :: ngrid(:)
    real(8), allocatable :: boxl(:,:)
    real(8), allocatable :: vpl(:,:)
    real(8), allocatable :: vpc(:,:)
    real(8), allocatable :: vpd(:)
    logical, allocatable :: mtpoint(:)
    integer, allocatable :: atom(:,:)
    integer, allocatable :: iv(:,:)
  end type rgrid

  private find_mt_points

contains

  !----------------------------------------------------------------------------
  subroutine print_rgrid(self)
    implicit none
    type(rgrid), intent(in) :: self
    integer :: ip, i, n
    write(*,*)
    write(*,*) "Box definition:"
    n = size(self%boxl,1)
    do i = 1, n
      write(*,'(3f12.4)') self%boxl(i,:)
    end do
    write(*,*) "Grid mesh:"
    write(*,*)  self%ngrid(:)
    write(*,*) "Real space grid:"
    write(*,*) " size: ", self%npt
    write(*,*) " ip  vpl  mtpoint  atom"
    do ip = 1, self%npt
      write(*,'(i4,3f12.4,L)') ip, self%vpl(:,ip), self%mtpoint(ip)
      if (self%mtpoint(ip)) then 
        write(*,'(2i4,4x,3i4)') self%atom(:,ip), self%iv(:,ip)
      end if
    end do
    write(*,*)
    return
  end subroutine

  !----------------------------------------------------------------------------
  subroutine delete_rgrid(self)
    implicit none
    type(rgrid), intent(inout) :: self
    if (allocated(self%ngrid))   deallocate(self%ngrid)
    if (allocated(self%boxl))    deallocate(self%boxl)
    if (allocated(self%vpl))     deallocate(self%vpl)
    if (allocated(self%vpc))     deallocate(self%vpc)
    if (allocated(self%vpd))     deallocate(self%vpd)
    if (allocated(self%mtpoint)) deallocate(self%mtpoint)
    if (allocated(self%atom))    deallocate(self%atom)
  end subroutine

  !----------------------------------------------------------------------------
  function gen_1d_rgrid(plot1d) result(self)
    implicit none
    ! input/output
    type(plot1d_type), intent(in) :: plot1d
    type(rgrid)                   :: self
    ! local
    integer :: iv, nv
    integer :: ip, np
    real(8), allocatable :: vvl(:,:)

    np = plot1d%path%steps
    nv = size(plot1d%path%pointarray)
    allocate(vvl(nv,3))
    do iv = 1, nv
      vvl(iv,:) = plot1d%path%pointarray(iv)%point%coord(:)
    end do
    ! path generator routine
    self = gen_1d(nv, vvl, np)
    deallocate(vvl)

    return
  end function

  !----------------------------------------------------------------------------
  function gen_1d(nv,vvl,np) result(self)
    implicit none
    ! input/output
    type(rgrid)         :: self
    integer, intent(in) :: nv         ! number of 1d path segments
    real(8), intent(in) :: vvl(nv,3)  ! segment coordinates
    integer, intent(in) :: np         ! number of grid points
    ! local
    integer :: iv
    integer :: ip, ip0, ip1, n
    real(8) :: t1, f, dt, vl(3), vc(3)
    real(8), allocatable :: seg(:), dv(:)

    ! initialization
    if (allocated(self%ngrid)) deallocate(self%ngrid)
    allocate(self%ngrid(1))
    self%ngrid(1) = np
    self%npt = np

    if (allocated(self%boxl)) deallocate(self%boxl)
    allocate(self%boxl(nv,3))
    self%boxl(:,:) = vvl(:,:)

    if (allocated(self%vpl)) deallocate(self%vpl)
    allocate(self%vpl(3,self%npt))
    if (allocated(self%vpc)) deallocate(self%vpc)
    allocate(self%vpc(3,self%npt))
    if (allocated(self%vpd)) deallocate(self%vpd)
    allocate(self%vpd(self%npt))

    ! find the total distance and the length of each segment
    allocate(dv(nv),seg(nv))
    dt = 0.d0
    do iv = 1, nv-1
      dv(iv) = dt
      vl(:) = vvl(iv+1,:)-vvl(iv,:)
      call r3mv(input%structure%crystal%basevect, vl, vc)
      seg(iv) = dsqrt(vc(1)**2+vc(2)**2+vc(3)**2)
      dt = dt+seg(iv)
    end do

    dv(nv) = dt
    if (dt < 1.d-8) Then
      do ip = 1, np
        self%vpl(:,ip) = vvl(1,:)
        self%vpd(ip) = 0.d0
      end do
    else
      do iv = 1, nv-1
        t1 = dble(np)*dv(iv)/dt
        ip0 = nint(t1)+1
        if (ip0 < 1) ip0 = 1
        t1 = dble(np)*dv(iv+1)/dt
        ip1 = nint(t1)
        if (ip1 > np) ip1 = np
        n = ip1-ip0
        if (n <= 0) n = 1
        do ip = ip0, ip1
          f = dble(ip-ip0)/dble(n)
          self%vpd(ip) = f*seg(iv)+dv(iv)
          self%vpl(:,ip) = vvl(iv,:)*(1.d0-f)+vvl(iv+1,:)*f
        end do
      end do
    end if
    deallocate(seg,dv)

    ! convert to cartesian coordinates
    do ip = 1, np
      call r3mv(input%structure%crystal%basevect, &
      &         self%vpl(:,ip), self%vpc(:,ip))
    end do

    ! find out which points belong to MT region
    call find_mt_points(self)

    return
  end function

  !----------------------------------------------------------------------------
  function gen_2d_rgrid(plot2d,iflag) result(self)
    implicit none
    ! input/output
    type(plot2d_type), intent(in) :: plot2d
    integer,           intent(in) :: iflag
    type(rgrid)                   :: self
    integer :: ngrid(2)
    real(8) :: boxl(3,3)

    ngrid(:)  = plot2d%parallelogram%grid(1:2)
    boxl(1,:) = plot2d%parallelogram%origin%coord(1:3)
    boxl(2,:) = plot2d%parallelogram%pointarray(1)%point%coord(1:3)-boxl(1,:)
    boxl(3,:) = plot2d%parallelogram%pointarray(2)%point%coord(1:3)-boxl(1,:)
    self = gen_2d(ngrid, boxl, iflag)

    return
  end function

  !----------------------------------------------------------------------------
  function gen_2d(ngrid,boxl,iflag) result(self)
    implicit none
    ! input/output
    type(rgrid)         :: self
    integer, intent(in) :: ngrid(2)
    real(8), intent(in) :: boxl(3,3)
    integer, intent(in) :: iflag
    ! local
    integer :: ip, ip1, ip2, iv(3)
    integer :: n1, n2, n0
    real(8) :: v1(3), v2(3), t1, t2

    if (allocated(self%ngrid)) deallocate(self%ngrid)
    allocate(self%ngrid(2))
    self%ngrid(:) = ngrid(:)

    if (allocated(self%boxl)) deallocate(self%boxl)
    allocate(self%boxl(3,3))
    self%boxl(:,:) = boxl(:,:)

    n1 = ngrid(1)
    n2 = ngrid(2)
    if (iflag>0) then
      ! periodic grid
      n0 = 1
      self%npt = n1*n2
    else
      ! general periodic grid (as employed by XCrySDen)
      n0 = 0
      self%npt = (n1+1)*(n2+1)
    end if

    if (allocated(self%vpl)) deallocate(self%vpl)
    allocate(self%vpl(3,self%npt))
    if (allocated(self%vpc)) deallocate(self%vpc)
    allocate(self%vpc(3,self%npt))

    ip = 0
    do ip2 = n0, n2
      t2 = dble(ip2-n0)/dble(n2)
      do ip1 = n0, n1
        ip = ip + 1
        t1 = dble(ip1-n0)/dble(n1)
        v1(:) = t1*boxl(2,:) &
        &     + t2*boxl(3,:) &
        &     + boxl(1,:)
        self%vpl(:,ip) = v1(:)
        call r3mv(input%structure%crystal%basevect, v1, v2)
        self%vpc(:,ip) = v2(:)
      end do
    end do

    ! find out which points belong to MT region
    call find_mt_points(self)

    return
  end function

  !----------------------------------------------------------------------------
  function gen_3d_rgrid(plot3d,iflag) result(self)
    implicit none
    ! input/output
    type(rgrid)                   :: self
    type(plot3d_type), intent(in) :: plot3d
    integer,           intent(in) :: iflag
    ! local
    integer :: ip, ip1, ip2, ip3, iv(3)
    integer :: n1, n2, n3, n0
    real(8) :: v1(3), v2(3), t1, t2, t3
    integer :: ngrid(3)
    real(8) :: boxl(4,3)

    ngrid(:)  = plot3d%box%grid(:)
    boxl(1,:) = plot3d%box%origin%coord(1:3)
    boxl(2,:) = plot3d%box%pointarray(1)%point%coord(1:3)-boxl(1,:)
    boxl(3,:) = plot3d%box%pointarray(2)%point%coord(1:3)-boxl(1,:)
    boxl(4,:) = plot3d%box%pointarray(3)%point%coord(1:3)-boxl(1,:)
    self = gen_3d(ngrid, boxl, iflag)

    return
  end function

  !----------------------------------------------------------------------------
  function gen_3d(ngrid,boxl,iflag) result(self)
    implicit none
    ! input/output
    type(rgrid) :: self
    integer, intent(in) :: ngrid(3)
    real(8), intent(in) :: boxl(4,3)
    integer, intent(in) :: iflag
    ! local
    integer :: ip, ip1, ip2, ip3, iv(3)
    integer :: n1, n2, n3, n0
    real(8) :: v1(3), v2(3), t1, t2, t3

    if (allocated(self%ngrid)) deallocate(self%ngrid)
    allocate(self%ngrid(3))
    self%ngrid(:) = ngrid(:)

    if (allocated(self%boxl)) deallocate(self%boxl)
    allocate(self%boxl(4,3))
    self%boxl(:,:) = boxl(:,:)

    n1 = ngrid(1)
    n2 = ngrid(2)
    n3 = ngrid(3)
    if (iflag>0) then
      ! periodic grid
      n0 = 1
      self%npt = n1*n2*n3
    else
      ! general periodic grid (as employed by XCrySDen)
      n0 = 0
      self%npt = (n1+1)*(n2+1)*(n3+1)
    end if

    if (allocated(self%vpl)) deallocate(self%vpl)
    allocate(self%vpl(3,self%npt))
    if (allocated(self%vpc)) deallocate(self%vpc)
    allocate(self%vpc(3,self%npt))

    ip = 0
    do ip3 = n0, n3
      t3 = dble(ip3-n0)/dble(n3)
      do ip2 = n0, n2
        t2 = dble(ip2-n0)/dble(n2)
        do ip1 = n0, n1
          t1 = dble(ip1-n0)/dble(n1)
          ip = ip+1
          v1(:) = t1*boxl(2,:) &
          &     + t2*boxl(3,:) &
          &     + t3*boxl(4,:) &
          &     + boxl(1,:)
          self%vpl(:,ip) = v1(:)
          ! convert point to Cartesian coordinates
          call r3frac(input%structure%epslat, v1, iv)
          call r3mv(input%structure%crystal%basevect, v1, v2)
          self%vpc(:,ip) = v2(:)
        end do
      end do
    end do

    ! find out which points belong to MT region
    call find_mt_points(self)

    return
  end function


  !----------------------------------------------------------------------------
  subroutine find_mt_points(self)
  	use modmain, only  : nspecies, natoms, atposc, rmt
    implicit none
    ! input/output
    type(rgrid),       intent(inout) :: self
    ! local variables
    integer :: ip, ia, is, i1, i2, i3, iv(3)
    integer :: n1, n2, n3
    real(8) :: t1, v0(3), v1(3), v2(3), v3(3), rmt2

    if (allocated(self%mtpoint)) deallocate(self%mtpoint)
    allocate(self%mtpoint(self%npt))
    self%mtpoint(:) = .false.

    if (allocated(self%atom)) deallocate(self%atom)
    allocate(self%atom(2,self%npt))
    self%atom(:,:) = 0

    if (allocated(self%iv)) deallocate(self%iv)
    allocate(self%iv(3,self%npt))
    self%iv(:,:) = 0    

    do ip = 1, self%npt
      v0(:) = self%vpl(:,ip)
      call r3frac(1.d-6, v0, iv)
      call r3mv(input%structure%crystal%basevect, v0, v1)
      ! check if point is in a muffin-tin
      do is = 1, nspecies
        rmt2 = rmt(is)**2
        do ia = 1, natoms(is)
          v2(:) = v1(:)-atposc(:,ia,is)
          ! apply back and forth translations
          do i3 = -1, 1
          do i2 = -1, 1
          do i1 = -1, 1
            v3(:) = v2(:)+ &
            &       dble(i1)*input%structure%crystal%basevect(:,1)+ &
            &       dble(i2)*input%structure%crystal%basevect(:,2)+ &
            &       dble(i3)*input%structure%crystal%basevect(:,3)
            t1 = v3(1)**2+v3(2)**2+v3(3)**2
            if (t1<rmt2) then
              ! belong to MT
              self%mtpoint(ip) = .true.
              self%atom(1,ip) = is
              self%atom(2,ip) = ia
              self%vpc(:,ip) = v3(:)
              self%iv(:,ip) = iv(:)-(/i1,i2,i3/)
              go to 10
            end if            
          end do
          end do
          end do
        end do ! ia
      end do ! is
      10 continue
    end do ! ip    

    return
  end subroutine

  !----------------------------------------------------------------------------
  subroutine calc_zdata_rgrid(self,ik,wfmt,wfir,zdata)
    use modmain
    implicit none

    type(rgrid), intent(in)  :: self
    integer,     intent(in)  :: ik
    complex(8),  intent(in)  :: wfmt(lmmaxapw,nrmtmax,natmtot)
    complex(8),  intent(in)  :: wfir(ngrtot)
    complex(8),  intent(out) :: zdata(*)
    ! local
    integer :: ip0, ip, is, ia, ias, ir, i, j, lm, ig, igp
    real(8) :: vl(3), vc(3), v(3), r
    complex(8), allocatable :: zylm(:)
    ! interpolation variables
    integer :: np2, ir0, iv(3), nx, ny, nz
    real(8) :: t1, t2, phs, phsav
    complex(8) :: zsum, z
    real(8), allocatable :: xa(:), ya(:), c(:)
    real(8), external :: polynom

    np2 = input%groundstate%nprad/2

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ip,is,ia,ias,vl,vc,zylm,r,ir,ir0,zsum,lm,j,i,xa,ya,t1,t2,c,v,igp)
#endif
    allocate(zylm(lmmaxapw))
    allocate(xa(input%groundstate%nprad))
    allocate(ya(input%groundstate%nprad))
    allocate(c(input%groundstate%nprad))
#ifdef USEOMP    
!$OMP DO
#endif
    do ip = 1, self%npt
      
      vl(:) = self%vpl(:,ip)
      vc(:) = self%vpc(:,ip)

    	if (self%mtpoint(ip)) then
        !----------------------------
        ! point belong to MT region
        !----------------------------
        ! cart. coordinates wrt atom (is,ia)
        is = self%atom(1,ip)
        ia = self%atom(2,ip)
        ias = idxas(ia,is)
        
        call ylm(vc,input%groundstate%lmaxapw,zylm)
        r = dsqrt(vc(1)**2+vc(2)**2+vc(3)**2)
        do ir = 1, nrmt(is)
          if (spr(ir,is)>=r) then
            if (ir<=np2) then
              ir0 = 1
            else if (ir>nrmt(is)-np2) then
              ir0 = nrmt(is)-input%groundstate%nprad+1
            else
              ir0 = ir-np2
            end if
            r = max(r,spr(1,is))
            zsum = 0.d0
            do lm = 1, lmmaxapw
            	do j = 1, input%groundstate%nprad
                i = ir0+j-1
                xa(j) = dble(wfmt(lm,i,ias))
                ya(j) = aimag(wfmt(lm,i,ias))
              end do
              ! interpolate radial part of wfmt
              t1 = polynom(0, input%groundstate%nprad, &
              &            spr(ir0,is), xa, c, r)
              t2 = polynom(0, input%groundstate%nprad, &
              &            spr(ir0,is), ya, c, r)
              zsum = zsum+cmplx(t1,t2,8)*zylm(lm)
            end do ! lm
            exit ! the loop over ir
          end if
        end do ! ir

        zdata(ip) = zsum

        v(:) = dble(self%iv(:,ip))
        t1 = twopi*(vkl(1,ik)*v(1)+ &
        &           vkl(2,ik)*v(2)+ &
        &           vkl(3,ik)*v(3))
        zdata(ip) = zsum*cmplx(dcos(t1),dsin(t1),8)

      else
        !----------------------------
        ! point belong to IS region
        !----------------------------
        zsum = 0.d0
        do igp = 1, ngk(1,ik)
          t1 = twopi*(vgkl(1,igp,1,ik)*vl(1)+ &
          &           vgkl(2,igp,1,ik)*vl(2)+ &
          &           vgkl(3,igp,1,ik)*vl(3))
          zsum = zsum + wfir(igp)*cmplx(dcos(t1),dsin(t1),8)
        end do
        zdata(ip) = zsum

      end if

    end do ! ip
#ifdef USEOMP
!$OMP END DO NOWAIT
#endif
    deallocate(zylm)
    deallocate(xa)
    deallocate(ya)
    deallocate(c)
#ifdef USEOMP
!$OMP END PARALLEL
#endif

    ! dephasing (lapw7)
    if (.false.) then
      phsav = 0.0d0
      do ip = 1, self%npt
        z = zdata(ip)
        if (abs(z)>1d-8) then
          phsav = phsav + dmod(datan2(dimag(z),dble(z))+pi,pi)
        end if
      end do
      phsav = phsav/dble(self%npt)
      phs   = dcmplx(cos(phsav),-sin(phsav))
      do ip = 1, self%npt
        zdata(ip) = zdata(ip)*phs
      end do
    end if

    ! t1 = dsqrt(0.5d0/omega)
    ! zdata(1:self%npt) = t1*zdata(1:self%npt)

    return
  end subroutine
  !----------------------------------------------------------------------------
  subroutine calc_zdata_rgrid_core(self,mu,corewfmt,zdata)
    use modmain
    use modxas
    implicit none

    type(rgrid), intent(in)  :: self
    integer,     intent(in)  :: mu
	real(8), intent(in) ::      corewfmt(spnrmax)
    complex(8),  intent(out) :: zdata(*)
    ! local
    integer :: ip0, ip, is, ia, ias, ir, i, j, lm, ig, igp ,lm1, lm2
    real(8) :: vl(3), vc(3), v(3), r
    complex(8), allocatable :: zylm(:)
    ! interpolation variables
    integer :: np2, ir0, iv(3), nx, ny, nz
    real(8) :: t1, t2, phs, phsav
    complex(8) :: zsum, z
    real(8), allocatable :: xa(:), ya(:), c(:)
    real(8), external :: polynom

    np2 = input%groundstate%nprad/2

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ip,is,ia,ias,vl,vc,zylm,r,ir,ir0,zsum,lm,j,i,xa,ya,t1,t2,c,v,igp)
#endif
    allocate(zylm(lmmaxapw))
    allocate(xa(input%groundstate%nprad))
    allocate(ya(input%groundstate%nprad))
    allocate(c(input%groundstate%nprad))
#ifdef USEOMP    
!$OMP DO
#endif
    do ip = 1, self%npt

      is = self%atom(1,ip)
      ia = self%atom(2,ip)
      ias = idxas(ia,is)
      vl(:) = self%vpl(:,ip)
      vc(:) = self%vpc(:,ip)

    	if (self%mtpoint(ip)) then
        !----------------------------
        ! point belong to MT region
        !----------------------------
        ! cart. coordinates wrt atom (is,ia)
        call ylm(vc,input%groundstate%lmaxapw,zylm)
        r = dsqrt(vc(1)**2+vc(2)**2+vc(3)**2)
        do ir = 1, nrmt(is)
          if (spr(ir,is)>=r) then
            if (ir<=np2) then
              ir0 = 1
            else if (ir>nrmt(is)-np2) then
              ir0 = nrmt(is)-input%groundstate%nprad+1
            else
              ir0 = ir-np2
            end if
            r = max(r,spr(1,is))
            lm1=idxlm(lxas,mj2ml(mj(mu),1))
            lm2=idxlm(lxas,mj2ml(mj(mu),2))
            zsum = 0.d0
            	do j = 1, input%groundstate%nprad
                i = ir0+j-1
                xa(j)=corewfmt(i)                
              end do
              ! interpolate radial part of wfmt
              t1 = polynom(0, input%groundstate%nprad, &
              &            spr(ir0,is), xa, c, r)
              zsum = zsum+t1*(1/(sqrt(2.0d0)))*(zylm(lm1)+zylm(lm2))
            exit ! the loop over ir
          end if
        end do ! ir
        write(*,*) 'ip=', ip
		
        zdata(ip) = zsum

        !v(:) = dble(self%iv(:,ip))
        !t1 = twopi*(vkl(1,ik)*v(1)+ &
        !&           vkl(2,ik)*v(2)+ &
        !&           vkl(3,ik)*v(3))
        !zdata(ip) = zsum*cmplx(dcos(t1),dsin(t1),8)

      else
        !----------------------------
        ! point belong to IS region
        !----------------------------
        zsum = 0.d0
        !do igp = 1, ngk(1,ik)
        !  t1 = twopi*(vgkl(1,igp,1,ik)*vl(1)+ &
        ! &           vgkl(2,igp,1,ik)*vl(2)+ &
        ! &           vgkl(3,igp,1,ik)*vl(3))
        ! zsum = zsum + wfir(igp)*cmplx(dcos(t1),dsin(t1),8)
        !end do
        zdata(ip) = zsum

      end if

    end do ! ip
#ifdef USEOMP
!$OMP END DO NOWAIT
#endif
    deallocate(zylm)
    deallocate(xa)
    deallocate(ya)
    deallocate(c)
#ifdef USEOMP
!$OMP END PARALLEL
#endif

    ! dephasing (lapw7)
    if (.false.) then
      phsav = 0.0d0
      do ip = 1, self%npt
        z = zdata(ip)
        if (abs(z)>1d-8) then
          phsav = phsav + dmod(datan2(dimag(z),dble(z))+pi,pi)
        end if
      end do
      phsav = phsav/dble(self%npt)
      phs   = dcmplx(cos(phsav),-sin(phsav))
      do ip = 1, self%npt
        zdata(ip) = zdata(ip)*phs
      end do
    end if

    ! t1 = dsqrt(0.5d0/omega)
    ! zdata(1:self%npt) = t1*zdata(1:self%npt)

    return
  end subroutine

end module
