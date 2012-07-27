!BOP
! 
! !ROUTINE: calcbradket
!
! !INTERFACE:
      subroutine calcbradket()

! !DESCRIPTION:
!
! This subroutine calculates the set of integrals:
!
! \begin{equation}
!\langle LM\tilde{\lambda} | \tilde{\lambda'}\rangle_a\equiv%
!\int^{R^a_{MT}}_{0}{\upsilon_{aLM}(r^a)%
!\tilde{u}^*_{\lambda}(r^a,E_{\lambda})\tilde{u}_{\lambda'}(r^a,E_{\lambda'})\left(r^a\right)^2dr^a}
!\end{equation}
!
! where $\tilde{u}$ can be $u$, $\dot{u}$ or $u^{lo}$
!
! !USES:
      
      use modinput
      use modmain
      use modgw

! !LOCAL VARIABLES:

      implicit none

      integer(4) :: ia
      integer(4) :: ias
      integer(4) :: ilo1
      integer(4) :: ilo2
      integer(4) :: io1
      integer(4) :: io2
      integer(4) :: ir
      integer(4) :: irm
      integer(4) :: is
      integer(4) :: ist1
      integer(4) :: ist2
      integer(4) :: l1
      integer(4) :: l2
      integer(4) :: lini1
      integer(4) :: lini2
      integer(4) :: ll1
      integer(4) :: ll2
      integer(4) :: nrwf

      real(8) :: fr(nrmtmax)
      real(8) :: gr(nrmtmax) 
      real(8) :: cf(3,nrmtmax)

! !DEFINED PARAMETERS:
      character(4) :: ftype(3)
      character(4) :: utype(3)
      data ftype /'core','apw ','lo  '/ 
      data utype /'    ','dot ','ddot'/
                                     
! 
! !INTRINSIC ROUTINES: 
!
      intrinsic dexp
      intrinsic iabs
      intrinsic isign
! 
! !REVISION HISTORY:
! 
! Created Dic 2003
! Last Modified: May. 2006
! Revisited 5.05.2011 by DIN 
!
!EOP
!BOC
!     Initializations
      nrwf=max(spnstmax,input%groundstate%lmaxapw,nlomax)
      allocate(bradketc(natmtot,maxnmix,spnstmax,0:nrwf,maxapword,3))
      allocate(bradketa(natmtot,maxnmix,0:input%groundstate%lmaxapw,&
     &  maxapword,0:nrwf,maxapword,3))
      allocate(bradketlo(natmtot,maxnmix,nlomax,0:nrwf,maxapword,3))
      bradketc=0.0d0
      bradketa=0.0d0
      bradketlo=0.0d0
      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          if(debug)write(701,100)ias,trim(spname(is))
!         Loop over radial mixed functions
          do irm=1,nmix(ias)
!---------------------------------------------------------------------------!
!             the left function of the product is a core wf.
!---------------------------------------------------------------------------!
            io1=1
            do ist1=1,ncore(is)
              l1=spl(ist1,is)
!
!             the right function of the product is a core wf.
!
              io2=1
              do ist2=1,ncore(is)
                l2=spl(ist2,is)
!               Check the triangular rule for L, l1 and l2
                if((iabs(l1-l2).le.bigl(ias,irm)).and.              &
               &   (l1+l2.ge.bigl(ias,irm)))then
                  do ir=1,nrmt(is)
                    fr(ir)=umix(ias,irm,ir)*rwfcr(ir,1,ist1,ias)*   &
                 &         rwfcr(ir,1,ist2,ias)*spr(ir,is)
                  enddo ! ir
!                 Calculate the integral:
                  call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
                  bradketc(ias,irm,ist1,ist2,io2,1)=gr(nrmt(is))
                  if(debug)write(701,1)irm,bigl(ias,irm),l1,ftype(1),utype(io1), &
                 &                     l2,ftype(1),utype(io2),                      &
                 &                     bradketc(ias,irm,ist1,ist2,io2,1)
                endif
              enddo ! ist2
!
!             the right function of the product is a valence wf.
!
              do l2=0,input%groundstate%lmaxapw
!               Check the triangular rule for L,l1 and l2
                if((iabs(l1-l2).le.bigl(ias,irm)).and.              &
               &   (l1+l2.ge.bigl(ias,irm)))then
                  do io2=1,apword(l2,is)
                    do ir=1,nrmt(is)
                      fr(ir)=umix(ias,irm,ir)*rwfcr(ir,1,ist1,ias)* &
                   &         apwfr(ir,1,io2,l2,ias)*spr(ir,is)
                    enddo ! ir
!                   Calculate the integral:
                    call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
                    bradketc(ias,irm,ist1,l2,io2,2)=gr(nrmt(is))
                    if(debug)write(701,1)irm,bigl(ias,irm),l1,ftype(1),      &
                   &                     utype(io1),l2,ftype(2),utype(io2),   &
                   &                     bradketc(ias,irm,ist1,l2,io2,2)
                  enddo ! io2
                endif
              enddo ! l2
!
!             the right function of the product is a local orbital.
!
              do ilo2=1,nlorb(is)
                l2=lorbl(ilo2,is)
                if((iabs(l1-l2).le.bigl(ias,irm)).and.              &
               &   (l1+l2.ge.bigl(ias,irm)))then
                    do io2=1,lorbord(ilo2,is)
                      do ir=1,nrmt(is)
                        fr(ir)=umix(ias,irm,ir)*rwfcr(ir,1,ist1,ias)*     &
                     &         lofr(ir,1,ilo2,ias)*spr(ir,is)
                      enddo ! ir
  !                   Calculate the integral:
                      call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
                      bradketc(ias,irm,ist1,ilo2,io2,3)=gr(nrmt(is))
                      if(debug)write(701,1)irm,bigl(ias,irm),l1,ftype(1),utype(io1),&
                     &                     l2,ftype(3),utype(io2),                  &
                     &                     bradketc(ias,irm,ist1,ilo2,io2,3)
                    end do ! io2
                endif
              enddo ! ilo2
            enddo ! ist1
!---------------------------------------------------------------------------!
!           the left function of the product is a valence wf.
!---------------------------------------------------------------------------!
            do l1=0,input%groundstate%lmaxapw
              do io1=1,apword(l1,is)
                
                lini1=ncore(is)+(io1-1)*input%groundstate%lmaxapw+1
                ll1=lini1+l1
!
!               the right function of the product is a core wf.
!
                io2=1
                do ist2=1,ncore(is)
                  l2=spl(ist2,is)
!                 Check the triangular rule for L, l1 and l2
                  if((iabs(l1-l2).le.bigl(ias,irm)).and.              &
                 &   (l1+l2.ge.bigl(ias,irm)))then
                    do ir=1,nrmt(is)
                      fr(ir)=umix(ias,irm,ir)*apwfr(ir,1,io1,l1,ias)* &
                   &         rwfcr(ir,1,ist2,ias)*spr(ir,is)
                    enddo ! ir
!                   Calculate the integral:
                    call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
                    bradketa(ias,irm,l1,io1,ist2,io2,1)=gr(nrmt(is))
                    if(debug)write(701,1)irm,bigl(ias,irm),l1,ftype(2),        &
                   &                     utype(io1),l2,ftype(1),utype(io2),     &
                   &                     bradketa(ias,irm,l1,io1,ist2,io2,1)
                  endif
                enddo ! ist2
!
!               the right function of the product is a valence wf.
!
                do l2=0,input%groundstate%lmaxapw
!                 Check the triangular rule for L,l1 and l2
                  if((iabs(l1-l2).le.bigl(ias,irm)).and.               &
                 &   (l1+l2.ge.bigl(ias,irm)))then
                    do io2=1,apword(l2,is)
                      do ir=1,nrmt(is)
                        fr(ir)=umix(ias,irm,ir)*apwfr(ir,1,io1,l1,ias)*&
                     &         apwfr(ir,1,io2,l2,ias)*spr(ir,is)
                      enddo ! ir
!                     Calculate the integral:
                      call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
                      bradketa(ias,irm,l1,io1,l2,io2,2)=gr(nrmt(is))
                      if(debug)write(701,1)irm,bigl(ias,irm),l1,ftype(2),         &
                     &                     utype(io1),l2,ftype(2),utype(io2),     &
                     &                     bradketa(ias,irm,l1,io1,l2,io2,2)
                    enddo ! io2
                  endif
                enddo ! l2
!
!               the right function of the product is a local orbital.
!
                do ilo2=1,nlorb(is)
                  l2=lorbl(ilo2,is)
                  if((iabs(l1-l2).le.bigl(ias,irm)).and.                &
                 &   (l1+l2.ge.bigl(ias,irm)))then
                    do io2=1,lorbord(ilo2,is)
                       do ir=1,nrmt(is)
                         fr(ir)=umix(ias,irm,ir)*apwfr(ir,1,io1,l1,ias)*  &
                      &         lofr(ir,1,ilo2,ias)*spr(ir,is)
                       enddo ! ir
!                      Calculate the integral:
                       call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
                       bradketa(ias,irm,l1,io1,ilo2,io2,3)=gr(nrmt(is))
                       if(debug)write(701,1)irm,bigl(ias,irm),l1,ftype(2),         &
                      &                     utype(io1),l2,ftype(3),utype(io2),      &
                      &                     bradketa(ias,irm,l1,io1,ilo2,io2,3)
                   end do !io2
                  endif
                enddo ! ilo2
              enddo ! io1 
            enddo ! l1 
!---------------------------------------------------------------------------!
!           The left function is a lo        
!---------------------------------------------------------------------------!
            do ilo1=1,nlorb(is)
              l1=lorbl(ilo1,is)
              ll1=ncore(is)+(input%groundstate%lmaxapw+1)*apwordmax+ilo1
              do io1=1,lorbord(ilo1,is)
!
!                the right function of the product is a core wf.
!
                 io2=1
                 do ist2=1,ncore(is)
                   l2=spl(ist2,is)
!                  Check the triangular rule for L,l1 and l2
                   if((iabs(l1-l2).le.bigl(ias,irm)).and.                &
                  &   (l1+l2.ge.bigl(ias,irm)))then
                     do ir=1,nrmt(is)
                       fr(ir)=umix(ias,irm,ir)*lofr(ir,1,ilo1,ias)* &
                    &         rwfcr(ir,1,ist2,ias)*spr(ir,is)
                     enddo ! ir
!                    Calculate the integral:
                     call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
                     bradketlo(ias,irm,ilo1,ist2,io2,1)=gr(nrmt(is))
                     if(debug)write(701,1)irm,bigl(ias,irm),l1,ftype(3),           &
                    &                     utype(io1),l2,ftype(1),utype(io2),       &
                    &                     bradketlo(ias,irm,ilo1,ist2,io2,1)
                   endif
                 enddo ! ist2
!
!                the right function of the product is a valencd wf.
!
                 do l2=0,input%groundstate%lmaxapw
!                  Check the triangular rule for L,l1 and l2
                   if((iabs(l1-l2).le.bigl(ias,irm)).and.                  &
                  &   (l1+l2.ge.bigl(ias,irm)))then
                     do io2=1,apword(l2,is)
                       lini2=ncore(is)+(io2-1)*input%groundstate%lmaxapw+1
                       ll2=lini2+l2
                       do ir=1,nrmt(is)
                         fr(ir)=umix(ias,irm,ir)*lofr(ir,1,ilo1,ias)*      &
                      &         apwfr(ir,1,io2,l2,ias)*spr(ir,is)
                       enddo ! ir
!                      Calculate the integral:
                       call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
                       bradketlo(ias,irm,ilo1,l2,io2,2)=gr(nrmt(is))
                       if(debug)write(701,1)irm,bigl(ias,irm),l1,ftype(3),           &
                      &                     utype(io1),l2,ftype(2),utype(io2),       &
                      &                     bradketlo(ias,irm,ilo1,l2,io2,2)
                     enddo ! io2
                   endif
                 enddo ! l2
!
!                the right function of the product is a local orbital.
!
                 do ilo2=1,nlorb(is)
                   l2=lorbl(ilo2,is)
                   ll2=ncore(is)+(input%groundstate%lmaxapw+1)*apwordmax+ilo2
                   if((iabs(l1-l2).le.bigl(ias,irm)).and.                  &
                  &   (l1+l2.ge.bigl(ias,irm)))then
                      do io2=1,lorbord(ilo2,is)
                        do ir=1,nrmt(is)
                          fr(ir)=umix(ias,irm,ir)*lofr(ir,1,ilo1,ias)*        &
                       &         lofr(ir,1,ilo2,ias)*spr(ir,is)
                        enddo ! ir
!                       Calculate the integral:
                        call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
                        bradketlo(ias,irm,ilo1,ilo2,io2,3)=gr(nrmt(is))
                        if(debug)write(701,1)irm,bigl(ias,irm),l1,ftype(3),             &
                       &                     utype(io1),l2,ftype(3),utype(io2),         &
                       &                     bradketlo(ias,irm,ilo1,ilo2,io2,3)
                      end do !io2
                   endif
                 enddo ! ilo2 
              end do !io1
            enddo ! ilo1
          enddo ! irm
        enddo ! ia
      enddo ! is
  100 format(/,5x,'Integrals <v_(NL)u_(l1)|u_(l2)> for atom ',I4,", ", &
     &       A,/,13x,'N',3x,'L',2x,'l1',1x,'u_',4x,'l2',1x,'u_',8x,    &
     &       '<v u|u>')
    1 format(10x,3i4,a4,1x,a4,i4,a4,1x,a4,1pe19.11)
      return
      end subroutine calcbradket
!EOC

