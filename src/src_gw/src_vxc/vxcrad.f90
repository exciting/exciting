
subroutine vxcrad
    use modinput
    use modmain, only : nrmtmax, nspecies, nrmt, spr, natoms, idxas, apword, &
    &                   idxlm, apwfr, vxcmt, nlorb, lorbl, lofr
    use mod_vxc, only : vxcraa, vxcrloa, vxcrlolo
    implicit none
    ! local variables
    integer :: is,ia,ias,nr,ir
    integer :: l1,l2,l3,m2,lm2
    integer :: ilo,ilo1,ilo2,io,io1,io2,lmmaxvr
    real(8) :: t1
    ! automatic arrays
    real(8) :: r2(nrmtmax), fr(nrmtmax), gr(nrmtmax), cf(3,nrmtmax)
    integer, allocatable :: lfromlm(:),mfromlm(:)

    lmmaxvr=(input%groundstate%lmaxvr+1)**2
    allocate (lfromlm(lmmaxvr))
    allocate (mfromlm(lmmaxvr))
    Do l2 = 0, input%groundstate%lmaxvr
      Do m2 = - l2, l2
        lm2 = idxlm (l2, m2)
        lfromlm(lm2)=l2
        mfromlm(lm2)=m2
      End Do
    End Do

    ! begin loops over atoms and species
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE) SHARED(lmmaxvr,lfromlm,mfromlm,natoms,vxcmt,vxcraa,vxcrloa,vxcrlolo,nspecies,nrmt,idxas,idxlm,apwfr,spr,lofr,apword,lorbl,input,nlorb) PRIVATE(l1,io1,l3,io2,l2,m2,lm2,ir,fr,cf,gr,nr,is,ia,ias,r2,t1)
#endif
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
#ifdef USEOMP
!$OMP DO
#endif
                  do lm2=1,lmmaxvr
                    l2=lfromlm(lm2)
                    m2=mfromlm(lm2)
!                  do l2 = 0, input%groundstate%lmaxvr
!                    do m2 = -l2, l2
!                      lm2 = idxlm(l2,m2)
                      do ir = 1, nr
                        t1 = apwfr(ir,1,io1,l1,ias)*apwfr(ir,1,io2,l3,ias)*r2(ir)
                        fr(ir) = t1*vxcmt(lm2,ir,ias)
                      end do
                      call fderiv(-1,nr,spr(:,is),fr,gr,cf)
                      vxcraa(io1,l1,io2,l3,lm2,ias)=gr(nr)
!                    end do
                  end do
#ifdef USEOMP
!$OMP END DO NOWAIT
#endif
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
!              do l2 = 0, input%groundstate%lmaxvr
#ifdef USEOMP
!$OMP DO
#endif
                  do lm2=1,lmmaxvr
                    l2=lfromlm(lm2)
                    m2=mfromlm(lm2)

!                do m2 = -l2, l2
!                  lm2 = idxlm(l2,m2)
                  do ir = 1, nr
                    t1 = lofr(ir,1,ilo,ias)*apwfr(ir,1,io,l3,ias)*r2(ir)
                    fr(ir) = t1*vxcmt(lm2,ir,ias)
                  end do
                  call fderiv(-1,nr,spr(:,is),fr,gr,cf)
                  vxcrloa(ilo,io,l3,lm2,ias) = gr(nr)
                end do
#ifdef USEOMP
!$OMP END DO NOWAIT
#endif
!              end do
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
!            do l2 = 0, input%groundstate%lmaxvr
#ifdef USEOMP
!$OMP DO
#endif
                  do lm2=1,lmmaxvr
                    l2=lfromlm(lm2)
                    m2=mfromlm(lm2)

!              do m2 = -l2, l2
!                lm2 = idxlm(l2,m2)
                do ir = 1, nr
                  t1 = lofr(ir,1,ilo1,ias)*lofr(ir,1,ilo2,ias)*r2(ir)
                  fr(ir) = t1*vxcmt(lm2,ir,ias)
                end do
                call fderiv(-1,nr,spr(:,is),fr,gr,cf)
                vxcrlolo(ilo1,ilo2,lm2,ias) = gr(nr)
              end do
#ifdef USEOMP
!$OMP END DO NOWAIT
#endif
!            end do
          end do
        end do
        ! end loops over atoms and species
      end do
    end do
#ifdef USEOMP
!$OMP END PARALLEL
#endif 
    deallocate(mfromlm,lfromlm)
    return
end subroutine
!EOC
