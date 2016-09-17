!BOP
! 
!!ROUTINE: calcbradket
!
!!INTERFACE:
!
subroutine calcbradket(ia,is)
!
!!DESCRIPTION:
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
!!USES:
    use modinput
    use modmain
    use modgw
      
!!INPUT VARIABLES:
    implicit none
    integer(4), intent(in) :: ia
    integer(4), intent(in) :: is
    
!!LOCAL VARIABLES:
    integer(4) :: ias
    integer(4) :: ilo1, ilo2
    integer(4) :: io1, io2
    integer(4) :: ir, irm
    integer(4) :: ist1, ist2
    integer(4) :: l1, l2
    real(8) :: fr(nrmtmax)
    real(8) :: gr(nrmtmax) 
    real(8) :: cf(3,nrmtmax)

!!DEFINED PARAMETERS:
    character(4) :: ftype(3)
    character(4) :: utype(3)
    data ftype /'core','apw ','lo  '/ 
    data utype /'    ','dot ','ddot'/
 
!!REVISION HISTORY:
! 
! Created Dic 2003
! Last Modified: May. 2006
! Revisited August 2012 by DIN 
!
!EOP
!BOC
      
    ias = idxas(ia,is)
    if (input%gw%debug) write(fdebug,100) ias, trim(spname(is))

    ! Loop over radial mixed functions
    do irm = 1, nmix(ias)

      if (input%gw%coreflag.ne."vab") then
        !------------------------------------------------
        ! the left function of the product is a core wf
        !------------------------------------------------
        io1 = 1
        do ist1 = 1, ncore(is)
          l1 = spl(ist1,is)

          !------------------------------------------------
          ! the right function of the product is a core wf
          !------------------------------------------------
          io2 = 1
          do ist2 = 1, ncore(is)
            l2 = spl(ist2,is)
            ! check the triangular rule for L, l1 and l2
            if ((iabs(l1-l2)<=bigl(irm,ias)).and.(bigl(irm,ias)<=(l1+l2))) then
              do ir = 1, nrmt(is)
                fr(ir) = umix(ir,irm,ias)*ucore(ir,1,ist1,ias)* &
                &        ucore(ir,1,ist2,ias)*spr(ir,is)
              end do ! ir
              ! calculate the integral:
              call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
              bradketc(1,irm,ist1,ist2,io2,ias) = gr(nrmt(is))
              ! debugging info
              if (input%gw%debug) &
              &  write(fdebug,1) irm, bigl(irm,ias), l1, &
              &                  ftype(1), utype(io1), l2, &
              &                  ftype(1), utype(io2), &
              &                  bradketc(1,irm,ist1,ist2,io2,ias)
            end if ! triangular rule
          end do ! ist2

          !----------------------------------------------------
          ! the right function of the product is a valence wf
          !----------------------------------------------------
          do l2 = 0, input%groundstate%lmaxapw
            ! check the triangular rule for L,l1 and l2
            if ((iabs(l1-l2)<=bigl(irm,ias)).and.(bigl(irm,ias)<=(l1+l2))) then
              do io2 = 1, apword(l2,is)
                do ir = 1, nrmt(is)
                  fr(ir) = umix(ir,irm,ias)*ucore(ir,1,ist1,ias)* &
                  &        apwfr(ir,1,io2,l2,ias)*spr(ir,is)
                end do ! ir
                ! calculate the integral
                call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
                bradketc(2,irm,ist1,l2,io2,ias) = gr(nrmt(is))
                ! debugging info
                if (input%gw%debug) &
                &  write(fdebug,1) irm, bigl(irm,ias), l1, &
                &                  ftype(1), utype(io1), l2, &
                &                  ftype(2), utype(io2), &
                &                  bradketc(2,irm,ist1,l2,io2,ias)
              end do ! io2
            end if ! triangular rule
          end do ! l2

          !-------------------------------------------------------
          ! the right function of the product is a local orbital
          !-------------------------------------------------------
          io2 = 1
          do ilo2 = 1, nlorb(is)
            l2 = lorbl(ilo2,is)
            if ((iabs(l1-l2)<=bigl(irm,ias)).and.(bigl(irm,ias)<=(l1+l2))) then
              do ir = 1, nrmt(is)
                fr(ir) = umix(ir,irm,ias)*ucore(ir,1,ist1,ias)* &
                &        lofr(ir,1,ilo2,ias)*spr(ir,is)
              end do ! ir
              ! calculate the integral
              call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
              bradketc(3,irm,ist1,ilo2,io2,ias) = gr(nrmt(is))
              ! debugging info
              if (input%gw%debug) &
              &  write(fdebug,1) irm, bigl(irm,ias), l1, &
              &                  ftype(1), utype(io1), l2, &
              &                  ftype(3), utype(io2), &
              &                  bradketc(3,irm,ist1,ilo2,io2,ias)
            end if ! triangular rule
          end do ! ilo2
              
        end do ! ist1
      end if ! left function is a core wf

      !-----------------------------------------------
      ! left function of the product is a valence wf
      !-----------------------------------------------
      do l1 = 0, input%groundstate%lmaxapw
        do io1 = 1, apword(l1,is)
                
          if (input%gw%coreflag.ne."vab") then
            !---------------------------------------------
            ! right function of the product is a core wf
            !---------------------------------------------
            io2 = 1
            do ist2 = 1, ncore(is)
              l2 = spl(ist2,is)
              ! check the triangular rule for L, l1 and l2
              if ((iabs(l1-l2)<=bigl(irm,ias)).and.(bigl(irm,ias)<=(l1+l2))) then
                do ir = 1, nrmt(is)
                  fr(ir) = umix(ir,irm,ias)*apwfr(ir,1,io1,l1,ias)* &
                  &        ucore(ir,1,ist2,ias)*spr(ir,is)
                end do ! ir
                ! calculate the integral
                call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
                bradketa(1,irm,l1,io1,ist2,io2,ias) = gr(nrmt(is))
                ! debugging info
                if (input%gw%debug) &
                &  write(fdebug,1) irm, bigl(irm,ias), l1, &
                &                  ftype(2), utype(io1), l2, &
                &                  ftype(1), utype(io2), &
                &                  bradketa(1,irm,l1,io1,ist2,io2,ias)
              end if ! triangular rule
            end do ! ist2
          end if

          !-----------------------------------------------
          ! right function of the product is a valence wf
          !-----------------------------------------------
          do l2 = 0, input%groundstate%lmaxapw
            ! check the triangular rule for L,l1 and l2
            if ((iabs(l1-l2)<=bigl(irm,ias)).and.(bigl(irm,ias)<=(l1+l2))) then
              do io2 = 1, apword(l2,is)
                do ir = 1, nrmt(is)
                  fr(ir) = umix(ir,irm,ias)*apwfr(ir,1,io1,l1,ias)* &
                  &        apwfr(ir,1,io2,l2,ias)*spr(ir,is)
                end do ! ir
                ! calculate the integral
                call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
                bradketa(2,irm,l1,io1,l2,io2,ias) = gr(nrmt(is))
                ! debugging info
                if (input%gw%debug) &
                &  write(fdebug,1) irm, bigl(irm,ias), l1, &
                &                  ftype(2), utype(io1), l2, &
                &                  ftype(2), utype(io2), &
                &                  bradketa(2,irm,l1,io1,l2,io2,ias)
              end do ! io2
            end if ! triangular rule
          end do ! l2

          !--------------------------------------------------
          ! right function of the product is a local orbital
          !--------------------------------------------------
          io2 = 1
          do ilo2 = 1, nlorb(is)
            l2 = lorbl(ilo2,is)
            if ((iabs(l1-l2)<=bigl(irm,ias)).and.(bigl(irm,ias)<=(l1+l2))) then
              do ir = 1, nrmt(is)
                fr(ir) = umix(ir,irm,ias)*apwfr(ir,1,io1,l1,ias)*  &
                &        lofr(ir,1,ilo2,ias)*spr(ir,is)
              end do ! ir
              ! calculate the integral
              call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
              bradketa(3,irm,l1,io1,ilo2,io2,ias) = gr(nrmt(is))
              ! debugging info
              if (input%gw%debug) &
              &  write(fdebug,1) irm, bigl(irm,ias), l1, &
              &                  ftype(2), utype(io1), l2, &
              &                  ftype(3), utype(io2), &
              &                  bradketa(3,irm,l1,io1,ilo2,io2,ias)
            end if ! triangular rule
          end do ! ilo2
        
        end do ! io1 
      end do ! l1 

      !---------------------------------------------------------------------------!
      ! left function is a lo        
      !---------------------------------------------------------------------------!
      io1 = 1
      do ilo1 = 1, nlorb(is)
        l1 = lorbl(ilo1,is)

        if (input%gw%coreflag.ne."vab") then
          !--------------------------------------------
          ! right function of the product is a core wf
          !--------------------------------------------
          io2 = 1
          do ist2 = 1, ncore(is)
            l2 = spl(ist2,is)
            ! check the triangular rule for L,l1 and l2
            if ((iabs(l1-l2)<=bigl(irm,ias)).and.(bigl(irm,ias)<=(l1+l2))) then
              do ir = 1, nrmt(is)
                fr(ir) = umix(ir,irm,ias)*lofr(ir,1,ilo1,ias)* &
                &        ucore(ir,1,ist2,ias)*spr(ir,is)
              end do ! ir
              ! calculate the integral:
              call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
              bradketlo(1,irm,ilo1,ist2,io2,ias) = gr(nrmt(is))
              ! debugging info
              if (input%gw%debug) &
              & write(fdebug,1) irm, bigl(irm,ias), l1, &
              &                  ftype(3), utype(io1), l2, &
              &                  ftype(1), utype(io2), &
              &                  bradketlo(1,irm,ilo1,ist2,io2,ias)
            end if ! triangular rule
          end do ! ist2
        end if

        !------------------------------------------------
        ! right function of the product is a valencd wf
        !------------------------------------------------
        do l2 = 0, input%groundstate%lmaxapw
          ! check the triangular rule for L,l1 and l2
          if((iabs(l1-l2)<=bigl(irm,ias)).and.(bigl(irm,ias)<=(l1+l2))) then
            do io2 = 1, apword(l2,is)
              do ir = 1, nrmt(is)
                fr(ir) = umix(ir,irm,ias)*lofr(ir,1,ilo1,ias)* &
                &        apwfr(ir,1,io2,l2,ias)*spr(ir,is)
              end do ! ir
              ! calculate the integral
              call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
              bradketlo(2,irm,ilo1,l2,io2,ias) = gr(nrmt(is))
              ! debugging info
              if (input%gw%debug) &
              &  write(fdebug,1) irm, bigl(irm,ias), l1, &
              &                  ftype(3), utype(io1), l2, &
              &                  ftype(2), utype(io2), &
              &                  bradketlo(2,irm,ilo1,l2,io2,ias)
            end do ! io2
          end if
        end do ! l2
      
        !---------------------------------------------------
        ! right function of the product is a local orbital
        !---------------------------------------------------
        io2 = 1
        do ilo2 = 1, nlorb(is)
          l2 = lorbl(ilo2,is)
          if ((iabs(l1-l2)<=bigl(irm,ias)).and.(bigl(irm,ias)<=(l1+l2))) then
            do ir = 1, nrmt(is)
              fr(ir) = umix(ir,irm,ias)*lofr(ir,1,ilo1,ias)* &
              &        lofr(ir,1,ilo2,ias)*spr(ir,is)
            end do ! ir
            ! calculate the integral
            call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
            bradketlo(3,irm,ilo1,ilo2,io2,ias) = gr(nrmt(is))
            ! debugging info
            if (input%gw%debug) &
            &  write(fdebug,1) irm, bigl(irm,ias), l1, &
            &                  ftype(3), utype(io1), l2, &
            &                  ftype(3), utype(io2), &
            &                  bradketlo(3,irm,ilo1,ilo2,io2,ias)
          end if ! triangular rule
        end do ! ilo2
      
      end do ! ilo1

    end do ! irm

    100 format(/,5x,'Integrals <v_(NL)u_(l1)|u_(l2)> for atom ',I4,", ", &
    &       A,/,13x,'N',3x,'L',2x,'l1',1x,'u_',4x,'l2',1x,'u_',8x,    &
    &       '<v u|u>')
    1 format(10x,3i4,a4,1x,a4,i4,a4,1x,a4,1pe19.11)
    
    return
end subroutine
!EOC

