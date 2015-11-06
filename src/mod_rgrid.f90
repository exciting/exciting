
module mod_rgrid
	use modinput
	implicit none

	type rgrid
		integer :: npt
		real(8), allocatable :: vpl(:,:)
		real(8), allocatable :: vpc(:,:)
		real(8), allocatable :: vpd(:)
		logical, allocatable :: mtpoint(:)
		integer, allocatable :: atom(:,:)
	end type rgrid

	private find_mt_points

contains

	!----------------------------------------------------------------------------
	subroutine print_rgrid(self)
		implicit none
		type(rgrid), intent(in) :: self
		integer :: ip
		write(*,*)
		write(*,*) "Real space grid:"
		write(*,*) " size: ", self%npt
		write(*,*) " ip  vpl  mtpoint  atom"
		do ip = 1, self%npt
			write(*,*) ip, self%vpl(:,ip), self%mtpoint(ip), self%atom(:,ip)
		end do
		return
	end subroutine

	!----------------------------------------------------------------------------
	subroutine delete_rgrid(self)
		implicit none
		type(rgrid), intent(inout) :: self
		if (allocated(self%vpl))     deallocate(self%vpl)
		if (allocated(self%vpc))     deallocate(self%vpc)
		if (allocated(self%vpd))     deallocate(self%vpd)
		if (allocated(self%mtpoint)) deallocate(self%mtpoint)
		if (allocated(self%atom))    deallocate(self%atom)
	end subroutine

	!----------------------------------------------------------------------------
  subroutine gen_1d_rgrid(self)
  	implicit none
  	! input/output
  	type(rgrid), intent(out) :: self
  	! local
  	integer :: nv, np, iv
  	integer :: ip, ip0, ip1, n
  	real(8) :: t1, f, dt, vl(3), vc(3)
  	real(8), allocatable :: vvl(:,:), seg(:), dv(:)

  	type(plot1d_type), pointer :: plotdef
  	plotdef => input%properties%wfplot%plot1d

		! number of 1d path segments
		nv = size(plotdef%path%pointarray)
		if (nv < 1) then
      write (*,*)
      write (*, '("Error(mod_rgrid::gen_1d_rgrid): nv < 1 : ", I8)') nv
      write (*,*)
      stop
    end if

    ! number of grid points
    np = plotdef%path%steps
    If (np < nv) then
      write (*,*)
      write (*, '("Error(connect): np < nv : ", 2I8)') np, nv
      write (*,*)
      stop
    end if

    ! initialization
    self%npt = np
    if (allocated(self%vpl)) deallocate(self%vpl)
  	allocate(self%vpl(3,self%npt))
		if (allocated(self%vpc)) deallocate(self%vpc)
		allocate(self%vpc(3,self%npt))
		if (allocated(self%vpd)) deallocate(self%vpd)
		allocate(self%vpd(self%npt))

		! find the total distance and the length of each segment
		allocate(vvl(3,nv))
	  do iv = 1, nv
      vvl(:,iv) = plotdef%path%pointarray(iv)%point%coord
  	end do

    allocate(dv(iv),seg(nv))
    dt = 0.d0
    do iv = 1, nv-1
     	dv(iv) = dt
     	vl(:) = vvl(:,iv+1)-vvl(:,iv)
     	call r3mv(input%structure%crystal%basevect, vl, vc)
     	seg(iv) = dsqrt(vc(1)**2+vc(2)**2+vc(3)**2)
     	dt = dt+seg(iv)
    end do

    dv(nv) = dt
    if (dt < 1.d-8) Then
      do ip = 1, np
        self%vpl(:,ip) = vvl(:,1)
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
          self%vpl(:,ip) = vvl(:,iv)*(1.d0-f)+vvl(:,iv+1)*f
        end do
      end do
    end if
    deallocate(vvl,seg,dv)

    ! convert to cartesian coordinates
    do ip = 1, np
      call r3mv(input%structure%crystal%basevect, &
      &         self%vpl(:,ip), self%vpc(:,ip))
    end do

    ! find out which points belong to MT region
    call find_mt_points(self)

    return
  end subroutine

  !----------------------------------------------------------------------------
  subroutine gen_2d_rgrid(self)
  	implicit none
  	! input/output
  	type(rgrid), intent(out) :: self
  	! local
  	integer :: ip, ip1, ip2, iv(3)
  	integer :: n1, n2
  	real(8) :: box(2,3), v1(3), v2(3), t1, t2

		type(plot2d_type), pointer :: plotdef
  	plotdef => input%properties%wfplot%plot2d

  	n1 = plotdef%parallelogram%grid(1)
  	n2 = plotdef%parallelogram%grid(2)
		self%npt = (n1+1)*(n2+1)

		if (allocated(self%vpl)) deallocate(self%vpl)
  	allocate(self%vpl(3,self%npt))
		if (allocated(self%vpc)) deallocate(self%vpc)
		allocate(self%vpc(3,self%npt))

		box(1,:) = plotdef%parallelogram%pointarray(1)%point%coord - &
    &          plotdef%parallelogram%origin%coord
    box(2,:) = plotdef%parallelogram%pointarray(2)%point%coord - &
    &          plotdef%parallelogram%origin%coord

    ! test whether box is reasonable ?
    
    ip = 0
    do ip2 = 0, n2
      t2 = dble(ip2)/dble(n2)
      do ip1 = 0, n1
        ip = ip + 1
        t1 = dble(ip1)/dble(n1)
      	v1(:) = t1*box(1,:)+ &
      	&       t2*box(2,:)+ &
      	&		    plotdef%parallelogram%origin%coord
      	self%vpl(:,ip) = v1(:)
      	call r3mv(input%structure%crystal%basevect, v1, v2)
      	self%vpc(:,ip) = v2(:)
      end do
    end do

    ! find out which points belong to MT region
    call find_mt_points(self)

  	return  	
  end subroutine

	!----------------------------------------------------------------------------
  subroutine gen_3d_rgrid(self)
  	implicit none
  	! input/output
  	type(rgrid), intent(out) :: self
  	! local
  	integer :: ip, ip1, ip2, ip3, iv(3)
  	integer :: n1, n2, n3
  	real(8) :: box(3,3), v1(3), v2(3), t1, t2, t3

  	type(plot3d_type), pointer :: plotdef
  	plotdef => input%properties%wfplot%plot3d

  	n1 = plotdef%box%grid(1)
  	n2 = plotdef%box%grid(2)
  	n3 = plotdef%box%grid(3)
  	self%npt = (n1+1)*(n2+1)*(n3+1)

		if (allocated(self%vpl)) deallocate(self%vpl)
  	allocate(self%vpl(3,self%npt))
		if (allocated(self%vpc)) deallocate(self%vpc)
		allocate(self%vpc(3,self%npt))

		box(1,:) = plotdef%box%pointarray(1)%point%coord-plotdef%box%origin%coord
    box(2,:) = plotdef%box%pointarray(2)%point%coord-plotdef%box%origin%coord
    box(3,:) = plotdef%box%pointarray(3)%point%coord-plotdef%box%origin%coord

  	ip = 0
    do ip3 = 0, n3
      t3 = dble(ip3)/dble(n3)
      do ip2 = 0, n2
        t2 = dble(ip2)/dble(n2)
        do ip1 = 0, n1
          t1 = dble(ip1)/dble(n1)
          ip = ip+1
          v1(:) = t1*box(1,:)+ &
          &       t2*box(2,:)+ &
          &       t3*box(3,:)+ &
          &       plotdef%box%origin%coord(:)
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
  end subroutine

  !----------------------------------------------------------------------------
  subroutine find_mt_points(self)
  	use modmain, only  : nspecies, natoms, atposc, rmt
  	implicit none
  	! input/output
  	type(rgrid),       intent(inout) :: self
  	! local variables
		integer :: ip, ia, is, i1, i2, i3, iv(3)
		real(8) :: t1, v1(3), v2(3), v3(3), rmt2

		if (allocated(self%mtpoint)) deallocate(self%mtpoint)
		allocate(self%mtpoint(self%npt))
		self%mtpoint(:) = .false.

		if (allocated(self%atom)) deallocate(self%atom)
		allocate(self%atom(2,self%npt))
		self%atom(:,:) = 0

		do ip = 1, self%npt
      v1(:) = self%vpc(:,ip)
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
	subroutine calc_zdata_rgrid(self,wfmt,wfir,zdata)
		use modmain, only : lmmaxvr, idxas, nrmt, spr, ngvec, vgc, &
		&                   natmtot, ngrtot, nrmtmax
		implicit none

	  type(rgrid), intent(in)  :: self
		complex(8),  intent(in)  :: wfmt(lmmaxvr,nrmtmax,natmtot)
    complex(8),  intent(in)  :: wfir(ngrtot)
    complex(8),  intent(out) :: zdata(*)
    ! local
    integer :: ip, is, ia, ias, ir, i, j, lm, ig, igp
    real(8) :: vl(3), vc(3), r
    complex(8), allocatable :: zylm(:)
	  ! interpolation variables
	  integer :: np2, ir0
	  real(8) :: t1, t2
	  complex(8) :: zsum
		real(8), allocatable :: xa(:), ya(:), c(:)
		real(8), external :: polynom

    allocate(zylm(lmmaxvr))
    allocate(xa(input%groundstate%nprad))
    allocate(ya(input%groundstate%nprad))
    allocate(c(input%groundstate%nprad))
    np2 = input%groundstate%nprad/2

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
    		call ylm(vc,input%groundstate%lmaxvr,zylm)
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
            do lm = 1, lmmaxvr
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

      else
				!----------------------------
    		! point belong to IS region
    		!----------------------------
    		zsum = 0.d0
        do ig = 1, ngvec
          t1 = vgc(1,ig)*vc(1)+ &
          &    vgc(2,ig)*vc(2)+ &
          &    vgc(3,ig)*vc(3)
          zsum = zsum + wfir(ig)*cmplx(dcos(t1),dsin(t1),8)
        end do
        zdata(ip) = zsum

      end if

    end do ! ip

    return
  end subroutine

end module