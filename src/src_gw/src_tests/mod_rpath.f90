
module mod_rpath

    use modinput
    implicit none

    type r_path
        integer :: atom1(2)
        integer :: atom2(2)
        integer :: nptot                    ! total number of points
        integer :: nptir                    ! number of points in interstitial region
        real(8) :: rvec(3)
        real(8) :: rlen
        real(8), allocatable :: r(:,:)      ! path coordinates
    end type r_path
    
    type(r_path) :: rpath
    
contains

    subroutine init_rpath(self,atom1,atom2)
        use modmain
        implicit none
        type(r_path), intent(out) :: self
        integer,     intent(in)  :: atom1, atom2
        ! local
        integer :: is, ia, is1, is2
        integer :: ir, jr, kr
        real(8) :: rlenir
    
        do is = 1, nspecies
        do ia = 1, natoms(is)
          if (idxas(ia,is)==atom1) then
            self%atom1(1) = is
            self%atom1(2) = ia
          endif  
        end do 
        end do
        
        if (atom1==atom2) then
          self%rvec(1:3) = input%structure%crystal%basevect(:,1)
          self%atom2(1) = self%atom1(1)
          self%atom2(2) = self%atom1(2)
        else  
          do is = 1, nspecies
          do ia = 1, natoms(is)
            if (idxas(ia,is)==atom2) then
              self%atom2(1) = is
              self%atom2(2) = ia
            end if  
          end do 
          end do  
          self%rvec(1:3) = atposc(1:3,self%atom2(2),self%atom2(1))- &
          &                atposc(1:3,self%atom1(2),self%atom1(1))
        end if  
        
        self%rlen = dsqrt(self%rvec(1)**2+self%rvec(2)**2+self%rvec(3)**2)
        
        ! Set real-space path from atom0 to atom1
        is1 = self%atom1(1)
        is2 = self%atom2(1)
        
        self%nptir = 100
        self%nptot = nrcmt(is1)+self%nptir+nrcmt(is2)-1
        
        allocate(self%r(self%nptot,1)) 
        do ir = 1, nrmt(is1)
          self%r(ir,1) = spr(ir,is1)
        end do
        rlenir = self%rlen-rmt(is1)-rmt(is2)
        do ir = 1, self%nptir-1
          jr = ir+nrmt(is1)
          self%r(jr,1) = rmt(is1)+dble(ir)*rlenir/dble(self%nptir)
        end do
        do ir = 1, nrmt(is2)
          jr = nrmt(is2)-ir+1
          kr = nrmt(is1)+self%nptir-1+ir
          self%r(kr,1) = self%rlen-spr(jr,is2)
        end do
        
        !write(77,*) "# self%nptot=", self%nptot
        !do ir = 1, self%nptot
        !  write(77,*) ir, self%r(ir,1)
        !end do
    
        return
    end subroutine

    !-----------------------------------------------------------
    subroutine init_radial_path(self,atom1)
        use modmain
        implicit none
        type(r_path), intent(out) :: self
        integer,      intent(in)  :: atom1
        ! local
        integer :: is, ia, is1, ia1, ir

        do is = 1, nspecies
        do ia = 1, natoms(is)
          if (idxas(ia,is)==atom1) then
            self%atom1(1) = is
            self%atom1(2) = ia
          endif  
        end do 
        end do

        is1 = self%atom1(1)
        ia1 = self%atom1(2)

        ! Direction
        self%rvec(1:3) = (/1.d0, 1.d0, 1.d0/)
        self%rlen = rmt(is1)
        
        ! Set a real-space radial path
        self%nptot = nrmt(is1)
        allocate(self%r(self%nptot,1))
        do ir = 1, nrmt(is1)
          self%r(ir,1) = spr(ir,is1)
        end do
        
        !write(77,*) "# self%nptot=", self%nptot
        !do ir = 1, self%nptot
        !  write(77,*) ir, self%r(ir)
        !end do
    
        return
    end subroutine

    !-----------------------------------------------------------
    subroutine init_azimuthal_path(self,atom1,ir0)
        use modmain
        implicit none
        type(r_path), intent(out) :: self
        integer,      intent(in)  :: atom1
        integer,      intent(in)  :: ir0

        ! local
        integer :: is, ia, is1, ir
        real(8) :: r, theta, phi

        do is = 1, nspecies
        do ia = 1, natoms(is)
          if (idxas(ia,is)==atom1) then
            self%atom1(1) = is
            self%atom1(2) = ia
          endif  
        end do 
        end do
        is1 = self%atom1(1)

        ! Spherical coordinates / Position of zero
        r = spr(ir0,is1)
        self%rlen = r
        phi = 0.d0

        ! Set a real-space path [-pi,pi] with a step 1 deg
        self%nptot = 361
        allocate(self%r(self%nptot,3))

        ! azimuthal path: [0, pi]
        do ir = 1, self%nptot
          theta = dble(ir-1)*twopi/dble(self%nptot)
          self%r(ir,1) = r*dsin(theta)*dcos(phi)
          self%r(ir,2) = r*dsin(theta)*dsin(phi)
          self%r(ir,3) = r*dcos(theta)
        end do

if (.false.) then
        write(85,*) "# self%nptot=", self%nptot
        do ir = 1, self%nptot
          theta = -pi+dble(ir-1)*twopi/dble(self%nptot)
          write(85,'(i4,4f16.6)') ir, theta*180.d0/pi, self%r(ir,:)
        end do
end if        
    
        return
    end subroutine


    !-----------------------------------------------------------
    subroutine calc_radial_wfmb(iq,npt,wfmb)
        use modmain
        use modgw
        implicit none

        integer, intent(in) :: iq
        integer, intent(in) :: npt
        complex(8), intent(out) :: wfmb(matsizmax,npt)
    
        integer :: is1, ia1, ias1, is2, ia2, ias2
        integer :: igq, jgq, ir, jr, kr, l, m, lm, imix, im, irm
        real(8) :: t1, ri(3), phs, phsat
        complex(8) :: eph
        complex(8), allocatable :: zylm(:)

        wfmb(:,:) = 0.d0
    
        ! calculate the MB functions
        is1 = rpath%atom1(1)
        ia1 = rpath%atom1(2)
        ias1 = idxas(ia1,is1)
    
        ! MT 1
        phsat = 2.0d0*pi*(kqset%vql(1,iq)*atposl(1,ia1,is1)+ &
                          kqset%vql(2,iq)*atposl(2,ia1,is1)+ &
                          kqset%vql(3,iq)*atposl(3,ia1,is1))
        eph = cmplx(cos(phsat),sin(phsat),8)

        allocate(zylm((maxbigl+1)*(maxbigl+1)))
        call ylm(rpath%rvec,maxbigl,zylm)
    
        im = 0
        do irm = 1, nmix(ias1)
          l = bigl(irm,ias1)  
          do m = -l, l
            im = im+1
            imix = locmixind(ias1,im)
            lm = idxlm(l,m)
            do ir = 1, rpath%nptot
              t1 = 1.0d0/spr(ir,is1)
              wfmb(imix,ir) = wfmb(imix,ir)+ &
              &               eph*cmplx(umix(ir,irm,ias1)*t1,0.0d0,8)*zylm(lm)
            end do
          end do
        end do
    
if (.false.) then
        if (iq==1) then
          do imix = 1, matsiz
            do ir = 1, npt
              write(87,'(i8,f18.6)') ir, dble(wfmb(imix,ir))
              write(88,'(i8,f18.6)') ir, imag(wfmb(imix,ir))
            end do
            write(87,*)
            write(88,*)
          end do
        end if
end if
    
        return   
    end subroutine

    !-----------------------------------------------------------
    subroutine calc_azimuthal_wfmb(iq,ir0,npt,wfmb)
        use modmain
        use modgw
        implicit none

        integer, intent(in) :: iq
        integer, intent(in) :: ir0
        integer, intent(in) :: npt
        complex(8), intent(out) :: wfmb(matsizmax,npt)
    
        integer :: is1, ia1, ias1, is2, ia2, ias2
        integer :: igq, jgq, ir, jr, kr, l, m, lm, imix, im, irm, lmmax
        real(8) :: t1, ri(3), phs, phsat
        complex(8) :: eph
        complex(8), allocatable :: zylm(:,:)

        wfmb(:,:) = 0.d0
    
        ! calculate the MB functions
        is1 = rpath%atom1(1)
        ia1 = rpath%atom1(2)
        ias1 = idxas(ia1,is1)
    
        ! MT 1
        phsat = 2.0d0*pi*(kqset%vql(1,iq)*atposl(1,ia1,is1)+ &
                          kqset%vql(2,iq)*atposl(2,ia1,is1)+ &
                          kqset%vql(3,iq)*atposl(3,ia1,is1))
        eph = cmplx(cos(phsat),sin(phsat),8)

        lmmax = (maxbigl+1)*(maxbigl+1)
        allocate(zylm(lmmax,rpath%nptot))
        do ir = 1, rpath%nptot
          call ylm(rpath%r(ir,1:3),maxbigl,zylm(:,ir))
        end do
    
        im = 0
        do irm = 1, nmix(ias1)
          l = bigl(irm,ias1)  
          do m = -l, l
            im = im+1
            imix = locmixind(ias1,im)
            lm = idxlm(l,m)
            do ir = 1, rpath%nptot
              t1 = 1.0d0/spr(ir0,is1)
              ! wfmb(imix,ir) = wfmb(imix,ir)+ &
              ! &               eph*cmplx(umix(ir0,irm,ias1)*t1,0.0d0,8)*zylm(lm,ir)
              wfmb(imix,ir) = eph*cmplx(umix(ir0,irm,ias1)*t1,0.0d0,8)*zylm(lm,ir)
            end do
          end do
        end do
        deallocate(zylm)
    
if (.false.) then
        if (iq==1) then
          do imix = 1, matsiz
            do ir = 1, npt
              write(87,'(i8,f18.6)') ir, dble(wfmb(imix,ir))
              write(88,'(i8,f18.6)') ir, imag(wfmb(imix,ir))
            end do
            write(87,*)
            write(88,*)
          end do
        end if
end if
    
        return   
    end subroutine      

end module
