!BOP
! !ROUTINE: symmbmt
!
! !INTERFACE:

subroutine genmbrotmat(iq,isym)

! !DESCRIPTION:
!
! The space group symmetry operations in MB representation
! 
! !USES:
      
      use modinput
      use modmain
      use modgw

! !LOCAL VARIABLES:

      implicit none
      integer(4), intent(in) :: iq
      integer(4), intent(in) :: isym
      
      integer(4) :: ia, is, ias
      integer(4) :: ja, js, jas
      integer(4) :: ieq, lspl
      integer(4) :: irm, jrm
      integer(4) :: l1, l2, m1, m2
      integer(4) :: imix,jmix   
      integer(4) :: im,jm   
      integer(4) :: iqp
      
      integer(4) :: ig, igq1, igq2
      integer(4) :: s(3,3)
      integer(4) :: g0(3), g1(3), g2(3), g3(3)

      real(8)    :: c(3,3)
      real(8)    :: t1
      real(8)    :: tstart, tend
      
      real(8)    :: v(3), v1(3)
      
      complex(8) :: zt1
      complex(8), allocatable :: tmat1(:,:), tmat2(:,:)
      complex(8), allocatable :: sgip(:,:)
      
! !EXTERNAL ROUTINES: 
      
      complex(8), external :: getdlmm
      
! !INTRINSIC ROUTINES: 


! !REVISION HISTORY:
!
! Created October 2011 by DIN
!      
!EOP
!BOC
      call cpu_time(tstart)
      if(tstart.lt.0.0d0)write(fgw,*)'warning, tstart < 0'
      
      if(allocated(rotmat))deallocate(rotmat)
      allocate(rotmat(matsiz,matsiz))
      rotmat(:,:)=zzero
      
!     index of the irreducible q-point in the full (non-reduced) grid      
      iqp=idikpq(indkpq(iq,1),1) 

      lspl=lsplsymc(isym)
      s(:,:)=symlat(:,:,lspl)
      c(:,:)=symlatc(:,:,lspl)
      
     ! G^{q}_{R}
      v(:)=matmul(vql(:,iq),dble(s))-vql(:,iqp)
      call r3frac(input%structure%epslat,v,g0)
     
!----------------------------------
!              MT part
!----------------------------------
!
!     Loop over atoms:
!
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias=idxas(ia,is)
!
!         Index of the symmetry equivalent atom: R[r_eq]=r_a
!          
          ieq=ieqatom(ia,is,isym)
!
!         Loop over mixed functions:
!
          imix=0
          do irm = 1, nmix(ias)
            l1=bigl(ias,irm)
            do m1 = -l1, l1
              imix=imix+1
              im=locmixind(ias,imix)
!
!             Loop over atoms:
!
              do js = 1, nspecies
                do ja = 1, natoms(js)
                  jas=idxas(ja,js)
!
!                 Loop over mixed functions:
!
                  jmix=0
                  do jrm = 1, nmix(jas)
                    l2=bigl(jas,jrm)
                    do m2 = -l2, l2
                      jmix=jmix+1
                      jm=locmixind(jas,jmix)

                      ieq=ieqatom(ia,is,isym)
                      if((ieq.eq.jas).and.(irm.eq.jrm).and.(l1.eq.l2))then
                        t1=twopi*dot_product(v(:),atposl(:,ia,is))
                        zt1=cmplx(cos(t1),sin(t1),8)
                        t1=-twopi*dot_product(vql(:,iq),vtlsymc(:,isym))
                        zt1=zt1*cmplx(cos(t1),sin(t1),8)
                        rotmat(jm,im)=zt1*getdlmm(c,l1,m2,m1)
                      endif

                    enddo ! m1
                  enddo ! irm
                enddo ! ia
              enddo ! is 

            enddo ! m2
          enddo ! jrm
        enddo ! ja
      enddo ! jas

!-----------------------------------------------------------------
!       Calculation of the matrix elements between two IPW's
!-----------------------------------------------------------------
      
      allocate(tmat1(ngq(iq),ngq(iq)))
      tmat1(:,:)=zzero
      
      do igq1 = 1, ngq(iqp)   ! loop over q0+G
        do igq2 = 1, ngq(iq)   ! loop over q1+G'
        
          g1(:)=ivg(:,igqig(igq1,iqp))                ! get G vector of q0+G
          g2(:)=ivg(:,igqig(igq2,iq))                 ! get G' vector of q1+G'
          g3(:)=matmul(g2(:),dble(s))-g1(:)+g0(:)
          ig=ivgig(g3(1),g3(2),g3(3))                 ! index of (RG-G'+G_{R}) vector
          if((ig.lt.1).or.(ig.gt.ngrtot)) then
            write(fgw,*) 'ERROR(getmbrotmat): Wrong ig!'
            stop 'genmbrotmat'
          end if
          
          v1(:)=vgql(:,igq2,iq)
          t1=-twopi*dot_product(v1(:),vtlsymc(:,isym))
          zt1=cmplx(cos(t1),sin(t1),8)
          
          tmat1(igq2,igq1)=zt1*ipwint(ig)
          
        enddo ! igq1
      enddo ! igq2

      allocate(sgip(ngq(iq),ngq(iq)))
      allocate(tmat2(ngq(iq),ngq(iq)))

      if(iq.ne.iqp)then
        ! save sgi for iq
        tmat2(:,:)=sgi(:,:)
        ! calculate sgi_p for iqp
        call diagsgi(iqp)
        sgip(:,:)=sgi(:,:)
        ! restore sgi for iq
        sgi(:,:)=tmat2(:,:)
      else
        sgip(:,:)=sgi(:,:)
      end if
      
      call zgemm('n','n',ngq(iq),ngq(iq),ngq(iq), &
     &           zone,tmat1,ngq(iq),sgip,ngq(iq),zzero,tmat2,ngq(iq))
     
      call zgemm('c','n',ngq(iq),ngq(iq),ngq(iq), &
     &           zone,sgi,ngq(iq),tmat2,ngq(iq),zzero,tmat1,ngq(iq))

      do imix=1,ngq(iq)
        do jmix=1,ngq(iq)
          rotmat(locmatsiz+jmix,locmatsiz+imix)=tmat1(imix,jmix)
        enddo
      enddo
     
      deallocate(tmat1,tmat2,sgip)
      
      call cpu_time(tend)
      if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
      call write_cputime(fgw,tend-tstart, 'GENMBROTMAT')
     
      return
end subroutine
!EOC
