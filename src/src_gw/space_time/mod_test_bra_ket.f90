module mod_braket_radial
    use modinput
    use modmain
    use modgw

    implicit none
    private
    public :: braket_radial

    contains

    subroutine braket_radial(ia,is)
        !
        !!DESCRIPTION:
        !
        ! This subroutine calculates the set of integrals:
        !
        !\begin{equation}
        !  \langle LM\tilde{\lambda} | \tilde{\lambda'}\rangle_a\equiv%
        !  \int^{R^a_{MT}}_{0}{\upsilon_{aLM}(r^a)%
        !  \tilde{u}^*_{\lambda}(r^a,E_{\lambda})\tilde{u}_{\lambda'}(r^a,E_{\lambda'})\left(r^a\right)^2dr^a}
        !\end{equation}
        !
        ! where $\tilde{u}$ can be $u$, $\dot{u}$ or $u^{lo}$
        !
        !!USES:
     
          
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
        integer(4) :: maxidx
        real(8) :: fr(nrmtmax)
        real(8) :: gr(nrmtmax) 
        real(8) :: cf(3,nrmtmax)
        
        real(8), allocatable :: brackets(:,:,:,:,:,:,:,:)
    
        !!DEFINED PARAMETERS:
        character(4) :: ftype(3)
        character(4) :: utype(3)
        data ftype /'core','apw ','lo  '/ 
        data utype /'    ','dot ','ddot'/
     
          
        ias = idxas(ia,is)
        
        ! maxidx = max([  input%groundstate%lmaxapw,   max(ncore(:nspecies)),       max(nlorb(:nspecies)) ])
        maxidx = max( input%groundstate%lmaxapw, maxval(ncore(:nspecies)), maxval(nlorb(:nspecies)) )
        
        allocate(brackets(natmtot, maxval(nmix(:nspecies)), 3, 3, maxapword, maxapword, maxidx, maxidx), source=0.0_dp)
    
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
                  ! bradketc(1,irm,ist1,ist2,io2,ias) = gr(nrmt(is))
                  brackets(ias, irm, 1, 1, 1, 1, ist1, ist2) = gr(nrmt(is))
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
                    ! bradketc(2,irm,ist1,l2,io2,ias) = gr(nrmt(is))
                    brackets(ias, irm, 1, 3, 1, io2, ist1, l2) = gr(nrmt(is))
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
                  ! bradketc(3,irm,ist1,ilo2,io2,ias) = gr(nrmt(is))
                  brackets(ias, irm, 1, 2, 1, 1, ist1, ilo2) = gr(nrmt(is))
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
                    ! bradketa(1,irm,l1,io1,ist2,io2,ias) = gr(nrmt(is))
                    brackets(ias, irm, 3, 1, io1, 1, l1, ist2) = gr(nrmt(is))
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
                    ! bradketa(2,irm,l1,io1,l2,io2,ias) = gr(nrmt(is))
                    brackets(ias, irm, 3, 3, io1, io2, l1, l2) = gr(nrmt(is))
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
                  ! bradketa(3,irm,l1,io1,ilo2,io2,ias) = gr(nrmt(is))
                  brackets(ias, irm, 3, 2, io1, 1, l1, ilo2) = gr(nrmt(is))
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
                  ! bradketlo(1,irm,ilo1,ist2,io2,ias) = gr(nrmt(is))
                  brackets(ias, irm, 2, 1, 1, 1, ilo1, ist2) = gr(nrmt(is))
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
                  ! bradketlo(2,irm,ilo1,l2,io2,ias) = gr(nrmt(is))
                  brackets(ias, irm, 2, 3, 1, io2, ilo1, l2) = gr(nrmt(is))
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
                ! bradketlo(3,irm,ilo1,ilo2,io2,ias) = gr(nrmt(is))
                brackets(ias, irm, 2, 2, 1, 1, ilo1, ilo2) = gr(nrmt(is))
              end if ! triangular rule
            end do ! ilo2
          
          end do ! ilo1
    
        end do ! irm
    end subroutine
    
end module mod_braket_radial
