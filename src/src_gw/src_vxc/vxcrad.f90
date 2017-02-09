
subroutine vxcrad
    use modinput
    use modmain, only : nrmtmax, nspecies, nrmt, spr, natoms, idxas, apword, &
    &                   idxlm, apwfr, vxcmt, nlorb, lorbl, lofr
    use mod_vxc, only : vxcraa, vxcrloa, vxcrlolo
    implicit none
    ! local variables
    integer :: is,ia,ias,nr,ir
    integer :: l1,l2,l3,m2,lm2
    integer :: ilo,ilo1,ilo2,io,io1,io2
    real(8) :: t1
    ! automatic arrays
    real(8) :: r2(nrmtmax), fr(nrmtmax), gr(nrmtmax), cf(3,nrmtmax)

    ! begin loops over atoms and species
    do is = 1, nspecies
      nr = nrmt(is)
      do ir = 1, nr
         r2(ir) = spr(ir,is)*spr(ir,is)
      end do
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        !---------------------------!
        !     APW-APW integrals     !
        !---------------------------!
        do l1 = 0, input%groundstate%lmaxmat
          do io1 = 1, apword(l1,is)
            do l3 = 0, input%groundstate%lmaxapw
              do io2 = 1, apword(l3,is)
                if (l1>=l3) then
                  do l2 = 0, input%groundstate%lmaxvr
                    do m2 = -l2, l2
                      lm2 = idxlm(l2,m2)
                      do ir = 1, nr
                        t1 = apwfr(ir,1,io1,l1,ias)*apwfr(ir,1,io2,l3,ias)*r2(ir)
                        fr(ir) = t1*vxcmt(lm2,ir,ias)
                      end do
                      call fderiv(-1,nr,spr(:,is),fr,gr,cf)
                      vxcraa(io1,l1,io2,l3,lm2,ias)=gr(nr)
                    end do
                  end do
                end if
              end do
            end do
          end do
        end do
        !--------------------------------------!
        !     local-orbital-APW integrals      !
        !--------------------------------------!
        do ilo = 1, nlorb(is)
          l1 = lorbl(ilo,is)
          do l3 = 0, input%groundstate%lmaxmat
            do io = 1, apword(l3,is)
              do l2 = 0, input%groundstate%lmaxvr
                do m2 = -l2, l2
                  lm2 = idxlm(l2,m2)
                  do ir = 1, nr
                    t1 = lofr(ir,1,ilo,ias)*apwfr(ir,1,io,l3,ias)*r2(ir)
                    fr(ir) = t1*vxcmt(lm2,ir,ias)
                  end do
                  call fderiv(-1,nr,spr(:,is),fr,gr,cf)
                  vxcrloa(ilo,io,l3,lm2,ias) = gr(nr)
                end do
              end do
            end do
          end do
        end do
        !-----------------------------------------------!
        !     local-orbital-local-orbital integrals     !
        !-----------------------------------------------!
        do ilo1 = 1, nlorb(is)
          l1 = lorbl(ilo1,is)
          do ilo2 = 1, nlorb(is)
            l3 = lorbl(ilo2,is)
            do l2 = 0, input%groundstate%lmaxvr
              do m2 = -l2, l2
                lm2 = idxlm(l2,m2)
                do ir = 1, nr
                  t1 = lofr(ir,1,ilo1,ias)*lofr(ir,1,ilo2,ias)*r2(ir)
                  fr(ir) = t1*vxcmt(lm2,ir,ias)
                end do
                call fderiv(-1,nr,spr(:,is),fr,gr,cf)
                vxcrlolo(ilo1,ilo2,lm2,ias) = gr(nr)
              end do
            end do
          end do
        end do
        ! end loops over atoms and species
      end do
    end do
    return
end subroutine
!EOC
