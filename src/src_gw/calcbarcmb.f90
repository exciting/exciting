!BOP
!
! !ROUTINE: calcbarcns
!
! !INTERFACE:
      subroutine calcbarcmb(iq)

! !DESCRIPTION:
!
!This subroutine calculates the matrix of the bare coulomb potential according to equations.
!
! !USES:

      use modmain
      use modgw 
      use modmpi, only: rank     
!
! !INPUT PARAMETERS: 

      implicit none

      integer(4), intent(in) :: iq ! index of the q-point

!
! !LOCAL VARIABLES:

      integer(4) :: ia
      integer(4) :: ias
      integer(4) :: is
      integer(4) :: iipw
      integer(4) :: ijrm
      integer(4) :: imix
      integer(4) :: ipw
      integer(4) :: irm
      integer(4) :: ja
      integer(4) :: jas
      integer(4) :: js
      integer(4) :: jmix
      integer(4) :: jrm       ! Indexes the mixed basis functions 
!                               of atom jat
      integer(4) :: kpmin     ! Starting value of the index for G vectors. 
!                               =2 if iq=1 (if q=0, discard G=0)
      integer(4) :: l1        ! Angular momentum quantum number of the 
!                               atomic function irm
      integer(4) :: l2        ! Angular momentum quantum number of the 
!                               atomic function jrm
      integer(4) :: lm1
      integer(4) :: lm12      ! joint index for the structure constant
!                               \Sigma(l1+l2,m3+m4)
      integer(4) :: m1        ! z-component angular momentum quantum 
!                               number of the atomic function irm
      integer(4) :: m2        ! z-component angular momentum quantum 
!                               number of the atomic function jrm
      integer(4) :: ylsize    ! Size of the array for spherical harmonics
      
      real(8) :: gpr          ! Scalar product G.r
      real(8) :: prefac       ! multipicative prefactor
      real(8) :: qg1len       ! |qg1|
      real(8) :: tg           ! value of \tilde{g}
      real(8) :: minu

      real(8) :: tstart, tend
      real(8) :: t1,t2

      real(8), dimension(3) :: gvec ! G-vector of the reciprocal 
!                                     lattice.
      real(8), dimension(3) :: gvecl ! G-vector of the reciprocal 
!                                     lattice.
      real(8), dimension(3) :: qg1  ! The sum q+G_1
      real(8), dimension(3) :: qvec ! q-point for which the bare 
!                                     coulomb potential is calculated.
      real(8), allocatable :: rtlij(:,:,:,:)
 
      complex(8) :: exev      ! exact eigenvalues of the coulomb
!                                matrix
      complex(8) :: expg      ! e^(iG.r)       
      complex(8) :: stc
      
      complex(8) :: y((maxbigl+1)*(maxbigl+1))! spherical harmonics 
     
      complex(8), allocatable :: sph(:,:)     ! spherical harmonics for
!                                               all G's
      complex(8), allocatable :: mat1(:,:),mat2(:,:)
      complex(8), allocatable :: vtemp(:,:)

!     for diagonalization subroutine
      real(8) :: vl, vu, abstol
      integer :: il, iu, neval, lwork, info, lrwork, liwork
      complex(8), allocatable :: work(:)
      real(8),    allocatable :: rwork(:)
      integer,    allocatable :: iwork(:), isuppz(:)
 
! !EXTERNAL ROUTINES: 
      
      real(8)  gettildeg
      external calcjlam
      external sigma
      external ylm
      external zgemm   
      real(8), external :: dlamch     

! !INTRINSIC ROUTINES: 

      intrinsic abs
      intrinsic conjg
      intrinsic sqrt

!
! !REVISION HISTORY:
! 
! Created 16th. March 2004 by RGA
! Last modified 31. March 2005 by RGA
! Revisited : May 2011 by DIN
!
!EOP
!BOC      
      call cpu_time(tstart)
      if(tstart.lt.0.0d0)write(fgw,*)'warning, tstart < 0'

      if(allocated(barc))deallocate(barc)
      allocate(barc(matsiz,matsiz))
      barc=zzero
      if(allocated(sqbarc))deallocate(sqbarc)
      allocate(sqbarc(matsiz,matsiz))
      sqbarc=zzero

!----------------------------------
!     Calculate barc matrix
!----------------------------------

!     Allocate local arrays
      allocate(vtemp(1:ngbarc(iq),1:locmatsiz))
      allocate(mat2(1:ngq(iq),1:locmatsiz))
      allocate(rtlij(natmtot,natmtot,maxnmix,maxnmix))

      prefac=16.0*pi*pi*sqrt(vi)

      ylsize=(maxbigl+1)*(maxbigl+1)
      allocate(sph(ylsize,ngbarc(iq)))
!
!     Calculate the matrix sigma for the structure constants
!
      call cpu_time(t1)
      call sigma(iq,4*(input%gw%MixBasis%lmaxmb+1))
      call cpu_time(t2)
      if (rank==0) call write_cputime(fgw,t2-t1,'    SIGMA')

!     the cartesian coordinates of the q-point
      qvec(1:3)=vqc(1:3,iq)
!
!     Set kpmin: Starting value of the index for G vectors. =2 if iq=1 (if q=0, discard G=0)
!
      kpmin=1
      if (Gamma) then
        kpmin=2
        call barcq0
      endif
!
!     Calculate all the products rtl*rtl
!
      rtlij=0.0d0
      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          do js=1,nspecies
            do ja=1,natoms(js)
              jas=idxas(ja,js)
              do irm=1,nmix(ias)
                do jrm=1,nmix(jas)
                  rtlij(ias,jas,irm,jrm)=rtl(ias,irm)*rtl(jas,jrm)
                enddo ! jrm  
              enddo ! irm  
            enddo ! ja
          enddo ! js  
        enddo ! ia  
      enddo ! is
!
!     Loop over atoms:
!
      imix=0
      do is = 1, nspecies
        do ia =1, natoms(is)
          ias=idxas(ia,is)

!         calculate the matrix elements jlam
          call calcjlam(iq,is,ias,mbl(ias))
!
!         Calculate Y_lm(q+G) for all G
!
          do iipw=1,ngbarc(iq)
            gvec(1:3)=vgc(1:3,igqigb(iipw,iq))
            qg1(1:3)=qvec(1:3)+gvec(1:3)
            call ylm(qg1,maxbigl,y)
            sph(1:ylsize,iipw)=y(1:ylsize)
          enddo
!
!         Loop over mixed functions:
!
          do irm = 1, nmix(ias)
            l1=bigl(ias,irm)
            do m1=-l1,l1
              imix = imix + 1
              
              jmix=0  
              do js = 1, nspecies
                do ja =1, natoms(js)
                  jas=idxas(ja,js)
                     
                  do jrm = 1, nmix(jas)
                    l2=bigl(jas,jrm)
                    do m2=-l2,l2
                      jmix=jmix+1
                      
                      if((.not.Gamma).or.(l1.ne.0).or.(l2.ne.0))then

                        ! old (RGA's implementation)
                        !tg=gettildeg(l1,l2,-m1,m2)
                        !lm12=(l1+l2)*(l2+l1+1)-m1+m2+1
                        !minu=(-1.0d0)**m1
                        !stc=tg*sgm(ias,jas,lm12)*minu
                        !barc(imix,jmix)=rtlij(ias,jas,irm,jrm)*stc

                        ! new
                        if(jas.ge.ias)then
                          tg=gettildeg(l1,l2,-m1,m2)
                          lm12=(l1+l2)*(l2+l1+1)+m2-m1+1
                          minu=(-1.0d0)**m1
                          stc=tg*sgm(ias,jas,lm12)
                        else
                          tg=gettildeg(l1,l2,m1,-m2)
                          lm12=(l1+l2)*(l2+l1+1)-m2+m1+1
                          minu=(-1.0d0)**(m2+l1+l2)
                          stc=tg*conjg(sgm(ias,jas,lm12))
                        endif
                        barc(imix,jmix)=rtlij(ias,jas,irm,jrm)*minu*stc
                        
                        ! Second term of Eq.(56)
                        if((m1.eq.m2).and.(l1.eq.l2).and.(ias.eq.jas))then
                          if(jrm.ge.irm)then
                            ijrm=irm+(jrm*(jrm-1))/2
                          else
                            ijrm=jrm+(irm*(irm-1))/2
                          endif  
                          barc(imix,jmix)=barc(imix,jmix)+ &
                            cmplx(rrint(ias,ijrm)*4.0d0*pi/dble(2*l1+1),0.0d0,8)
                        end if
                      
                      end if
                      
                    enddo ! m2
                  enddo ! jrm
                enddo ! ja
              enddo ! js 
!
!             Calculation of the matrix element between an atomic
!             mixed function and an IPW
!
             if (Gamma) vtemp(1,imix)=zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ipw,gvecl,gpr,expg,qg1,qg1len,lm1)
!$OMP DO
#endif
             do ipw = kpmin, ngbarc(iq)

               gvecl(1:3)=dble(ivg(1:3,igqigb(ipw,iq)))
               gpr=gvecl(1)*atposl(1,ia,is)+ &
              &    gvecl(2)*atposl(2,ia,is)+ &
              &    gvecl(3)*atposl(3,ia,is)
               expg=cmplx(cos(2.0d0*pi*gpr),-sin(2.0d0*pi*gpr),8)

               qg1(1:3) = qvec(1:3) + vgc(1:3,igqigb(ipw,iq))
               qg1len = qg1(1)*qg1(1) + qg1(2)*qg1(2) + qg1(3)*qg1(3)
               
               lm1=l1*(l1+1)+m1+1
               vtemp(ipw,imix)=cmplx(prefac*jlam(irm,ipw)/qg1len,0.0d0,8)*  &
              &    sph(lm1,ipw)*((-zi)**l1)*expg
              
              enddo ! ipw
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            enddo ! m1
          enddo ! irm
        enddo ! ieq
      enddo ! iat
      
      call zgemm('n','n',ngq(iq),locmatsiz,ngbarc(iq),zone,mpwipw,ngq(iq),   &
     &    vtemp,ngbarc(iq),zzero,mat2,ngq(iq)) 

      do iipw=1,ngq(iq)
        do imix=1,locmatsiz
          barc(locmatsiz+iipw,imix)=mat2(iipw,imix)
          barc(imix,locmatsiz+iipw)=conjg(mat2(iipw,imix))
        enddo
      enddo    
      deallocate(mat2,jlam,rtlij,sph,vtemp)
      deallocate(sgm)
      
      allocate(mat1(1:ngq(iq),1:ngbarc(iq)))
      allocate(mat2(1:ngq(iq),1:ngq(iq)))
!
!     Calculation of the matrix elements between two IPW's
!
      mat1(:,1)=zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ipw,qg1,qg1len,exev,iipw)
!$OMP DO
#endif
      
      do ipw = kpmin, ngbarc(iq)
        qg1(1:3) = qvec(1:3) + vgc(1:3,igqigb(ipw,iq))
        qg1len = qg1(1)*qg1(1) + qg1(2)*qg1(2) + qg1(3)*qg1(3)
        exev = cmplx(4.0d0*pi/qg1len,0.0d0,8)
        do iipw=1,ngq(iq)
          mat1(iipw,ipw)=mpwipw(iipw,ipw)*exev
        enddo  
      enddo
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
      
      call zgemm('n','c',ngq(iq),ngq(iq),ngbarc(iq),zone,mat1,ngq(iq), &
     &    mpwipw,ngq(iq),zzero,mat2,ngq(iq)) 
      
      do ipw=1,ngq(iq)
        do iipw=1,ngq(iq)
          barc(locmatsiz+ipw,locmatsiz+iipw)=mat2(ipw,iipw)
        enddo  
      enddo
      deallocate(mat1)
      deallocate(mat2)
      
!----------------------------------
!     Diagonalize the bare coulomb matrix
!----------------------------------

      if(allocated(vmat))deallocate(vmat)
      allocate(vmat(matsiz,matsiz))
      
      if(allocated(barcev))deallocate(barcev)
      allocate(barcev(matsiz))

if (.false.) then      
      vmat(1:matsiz,1:matsiz) = barc(1:matsiz,1:matsiz)
      deallocate(barc)  
      lwork = 2*matsiz
      allocate(work(lwork),rwork(3*matsiz))
      call zheev( 'v','u',matsiz,vmat,matsiz, &
      &           barcev,work,lwork,rwork,info)
      call errmsg(info.ne.0,'CALCBARCMB',"Fail to diag. barc by zheev !!!")
      deallocate(work,rwork)
else
      lrwork = -1
      liwork = -1
      lwork = -1
      iu = matsiz
      abstol = 2.d0*dlamch('S')
      allocate(work(1),rwork(1),iwork(1))
      allocate(isuppz(2*matsiz))
      call zheevr('V', 'A', 'U', matsiz, barc, matsiz, vl, vu, il, iu, &
      &           abstol, neval, barcev, vmat, matsiz, isuppz, work, lwork, rwork, &
      &           lrwork, iwork, liwork, info)
      call errmsg(info.ne.0,'CALCBARCMB',"Fail to diag. barc by zheevr !!!")
      lrwork=int(rwork(1))
      liwork=int(iwork(1))
      lwork=int(work(1))
      ! write(*,*) lrwork,liwork,lwork
      deallocate(work,rwork,iwork)
      allocate(work(lwork),rwork(lrwork),iwork(liwork))
      call zheevr('V', 'A', 'U', matsiz, barc, matsiz, vl, vu, il, iu, &
      &           abstol, neval, barcev, vmat, matsiz, isuppz, work, lwork, rwork, &
      &           lrwork, iwork, liwork, info)
      call errmsg(info.ne.0,'CALCBARCMB',"Fail to diag. barc by zheevr !!!")
      deallocate(work,rwork,iwork,isuppz)
      deallocate(barc)
end if
     
      if(debug)then
        write(55,*) "### barc eigenvalues ###"
        do imix=1,matsiz
          write(55,'(i5,e16.6)') imix, barcev(imix) 
        enddo 
      endif !debug

      call cpu_time(tend)
      if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
      if (rank==0) call write_cputime(fgw,tend-tstart, 'CALCBARCMB')

      end subroutine calcbarcmb
!EOC      
