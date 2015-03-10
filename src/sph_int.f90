Subroutine sph_int
!single center integration

      Use modmain, Only: lmmaxvr, rhoir, rhomt
      Use modinput
      Implicit None
      Integer :: lmax, nsph, nr, ntot
      Integer :: ir, icount, jcount
      Integer :: idxmin, idxmax, idxmodulo
      Real(8) :: pi, rm, r, rmin, G, wr, f, I
      Real(8), allocatable :: v_unitsph(:,:), vpc(:,:), vpl(:,:)
      Real(8) :: vorigin(3)
      Real(8), allocatable :: wsph(:), x(:), fp(:)
      Real(8) :: inv_basevect(3,3)

!---------------------------------------------
! Read density from file
!---------------------------------------------
      Call init0
! read density from file
      Call readstate

      pi=atan(1.)*4
      vorigin=(/ 0, 0, 0 /)
      rm=1.4
      nsph=230 !possible numbers are: 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3740, 3890, 4334, 4802, 5294, 5810
      nr=30
      ntot=nsph*nr
      allocate(v_unitsph(3,nsph), vpc(3,ntot), vpl(3,ntot))
      allocate(wsph(nsph), x(nr), fp(ntot))
      lmax=input%groundstate%lmaxvr
      call leblaik(nsph, v_unitsph(:,:),wsph(:))
      wsph(:)=wsph(:)*4*pi
      do ir=1,nr
         x(ir)=cos(dble(ir)/(nr+1)*pi)
         r=rm*(1+x(ir))/(1-x(ir))
         idxmin=(ir-1)*nsph + 1
         idxmax=idxmin + nsph -1
         vpc(:,idxmin:idxmax) = r*v_unitsph(:,:)
      end do
      rmin=rm*(1+x(nr))/(1-x(nr))
      do jcount=1,3
         if (abs(vorigin(jcount)) .lt. epsilon(rmin)*rmin/100) then  
            do icount=1,ntot
               vpc(jcount,icount) = vpc(jcount,icount) + vorigin(jcount)
            enddo
         endif
      enddo

      call r3minv (input%structure%crystal%basevect, inv_basevect)
      vpl(:,:)=matmul(inv_basevect(:,:),vpc(:,:))
      call rfarray(lmax,lmmaxvr,rhomt,rhoir,ntot,vpl,fp)
      do icount=1,ntot
         if      ( (vpl(1,icount) .gt. 0.5) .or. (vpl(1,icount) .lt. -0.5) & 
              .or. (vpl(2,icount) .gt. 0.5) .or. (vpl(2,icount) .lt. -0.5) &
              .or. (vpl(3,icount) .gt. 0.5) .or. (vpl(3,icount) .lt. -0.5) ) then
           fp(icount) = 0.
         endif
      enddo

      I=0.
      do ir=1,nr
         idxmin=(ir-1)*nsph + 1
         idxmax=idxmin + nsph -1
         idxmodulo=0
         G=0.
         do icount = idxmin, idxmax
            idxmodulo = idxmodulo + 1
            G = G + wsph(idxmodulo)*fp(icount)
         enddo
         wr=pi/(nr+1)*(sin(dble(ir)/(nr+1)*pi))**2
         f=((1+x(ir))/(1-x(ir))**3)**(1.5)*G
         r=rm*(1+x(ir))/(1-x(ir))
!         if (r .le. 5) then
            I = I + wr*f
!         endif
      enddo
I = I*2*rm**3
write(*,*)'I = ', I
End Subroutine

