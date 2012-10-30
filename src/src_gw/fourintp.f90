!BOP
!
! !ROUTINE: fourintp
!
! !INTERFACE:
      subroutine fourintp(nb,nk1,kvecs1,f1,nk2,kvecs2,f2)
      
! !DESCRIPTION:
!
! This subroutine interpolate function $f1_n( k )$ defined on the kmesh 1, to kmesh 2
! using 3D Smooth Fourier transform according to 
! PRB 38, 2721 (1988).
!
! !USES:
       
      use modinput
      use modmain,   only: nsymcrys,pi,zzero,zone
      use modgw,     only: avec
      use fouri,     only: nrr,nst,rindex,rst,setrindex_done


! !LOCAL VARIABLES:

      implicit none
      integer(4),intent(in) :: nb,nk1,nk2
      real(8),   intent(in) :: kvecs1(3,nk1),kvecs2(3,nk2)
      complex(8),intent(in) :: f1(1:nb,nk1)
      complex(8),intent(out):: f2(1:nb,nk2) 
      
      integer(4) :: i
      integer(4) :: ib
      integer(4) :: ik
      integer(4) :: ir
      integer(4) :: info
      integer(4) :: ist
      integer(4) :: jk
      integer(4), allocatable  :: ipiv(:)
      
      real(8) :: den
      real(8) :: pref
      real(8) :: kdotr
      real(8) :: rmin,rlen,x2,x6,c1,c2
      real(8), dimension(3) :: r,rvec,kvec
      real(8), allocatable  :: rho(:)
      
      complex(8) :: expkr
      
      complex(8), allocatable :: dele(:,:)
      complex(8), allocatable :: h(:,:)
      complex(8), allocatable :: coef(:,:)
      complex(8), allocatable :: smat1(:,:),smat2(:,:)
      complex(8), allocatable :: sm2(:,:)

! !DEFINED PARAMETERS:
 
! !REVISION HISTORY:
!
! Created 02.03.06 by RGA
! Revisited July 2011 by DIN
!
!EOP
!BOC
      if (input%gw%debug) call linmsg(6,'-',' Fourier interpolation: Start')

!     shortcut for basis vectors 
      avec(:,1)=input%structure%crystal%basevect(:,1)
      avec(:,2)=input%structure%crystal%basevect(:,2)
      avec(:,3)=input%structure%crystal%basevect(:,3)
!
!     Set rindex to prepare for Fourier interpolation
!
      if(.not.setrindex_done) then  
        call setrindex
        setrindex_done = .true.
      endif 

      den=dble(nsymcrys)

      c1=0.25d0
      c2=0.25d0
      allocate(smat1(nst,nk1),     &
     &         smat2(nst,nk2),     &
     &         rho(nst),           &
     &         coef(nb,nst),       &
     &         ipiv(1:nk1-1),      &
     &         sm2(1:nst,1:nk1-1), &
     &         h(1:nk1-1,1:nk1-1), &
     &         dele(1:nb,1:nk1-1))
!
!     Calculate the star expansion function at each irreducible k-point
!
      smat1(1:nst,1:nk1)=zzero
      do ik=1,nk1
        kvec(1:3)=kvecs1(1:3,ik)
        do ir=2,nrr
          ist=rst(1,ir)
          pref=dble(rst(2,ir))
          r(1:3)=dble(rindex(1:3,ir))
          kdotr=2.0d0*pi*(r(1)*kvec(1)+r(2)*kvec(2)+r(3)*kvec(3)) 
          expkr=cmplx(cos(kdotr),sin(kdotr),8)
          smat1(ist,ik)=smat1(ist,ik)+pref*expkr/den
        enddo
      enddo 
!
!     Calculate the curvature function (rho) for each star
!
      rho(1:nst)=0.0d0
      ist=1
      do ir=2,nrr
        if(rst(1,ir).ne.ist)then
          ist=rst(1,ir)
          r(1:3)=dble(rindex(1:3,ir))
          do i=1,3
             rvec(i)=r(1)*avec(i,1)+ &
            &        r(2)*avec(i,2)+ &
            &        r(3)*avec(i,3)   
          end do
          rlen=rvec(1)*rvec(1)+rvec(2)*rvec(2)+rvec(3)*rvec(3)
          if(ist.eq.2) rmin=rlen
          x2=rlen/rmin
          x6=x2*x2*x2
          rho(ist)=(1.0d0-c1*x2)*(1.0d0-c1*x2)+c2*x6
        endif
      enddo
!
!      Set sm2(k)=smat(k)-smat(k_nkp) and dele
!
      do ik=1,nk1-1
        do ist=2,nst
          sm2(ist,ik)=smat1(ist,ik)-smat1(ist,nk1)
        enddo
        do ib=1,nb
          dele(ib,ik)=f1(ib,ik)-f1(ib,nk1)
        enddo
      enddo
!
!     Calculate the matrix H      
!
      h(1:nk1-1,1:nk1-1)=zzero
      do ik=1,nk1-1
        do jk=1,nk1-1
          do ist=2,nst
            h(ik,jk)=h(ik,jk)+sm2(ist,ik)*conjg(sm2(ist,jk))/rho(ist)
          enddo
        enddo
      enddo
!
!     Solve the Linear equations for the Lagrange multipliers
!
      call zgetrf(nk1-1,nk1-1,h,nk1-1,ipiv,info)
      call errmsg(info.ne.0,"fourintp","error when calling zgetrf")

      call zgetrs('n',nk1-1,nk1-1,h,nk1-1,ipiv,dele,nb,info)
      call errmsg(info.ne.0,"fourintp","error when calling zgetrs")
!
!     Calculate the coefficients of the Star expansion
!
      coef(1:nb,1)=f1(1:nb,nk1)
      do ist=2,nst
        coef(1:nb,ist)=zzero
        do ik=1,nk1-1
          coef(1:nb,ist)=coef(1:nb,ist)+dele(1:nb,ik)*conjg(sm2(ist,ik))
        enddo
        coef(1:nb,ist)=coef(1:nb,ist)/rho(ist)
        coef(1:nb,1)=coef(1:nb,1)-coef(1:nb,ist)*smat1(ist,nk1)
      enddo
!
!     Calculate the interpolated function on the new mesh 
!
      smat2(1:nst,1:nk2)=zzero
      do ik=1,nk2
        kvec(1:3)=kvecs2(1:3,ik)
        do ir=1,nrr
          ist=rst(1,ir)
          pref=dble(rst(2,ir))
          r(1:3)=dble(rindex(1:3,ir))
          kdotr=2.0d0*pi*sum(r*kvec) 
          expkr=cmplx(cos(kdotr),-sin(kdotr),8)
          smat2(ist,ik)=smat2(ist,ik)+pref*expkr/den
        enddo 
      enddo 
      call zgemm('n','n',nb,nk2,nst,zone,coef,nb,smat2,nst,zzero,f2,nb)
      
      if (input%gw%debug) call linmsg(6,'-','Fourier interpolation: Done')

      deallocate(smat1,smat2,coef,rho,ipiv,sm2,h,dele)
      
      return
      end subroutine fourintp
!EOC   
      
      
