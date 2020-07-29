
subroutine vmat_ldapu(is, ia, ngk, apwalm, evecfv, vnn)
    use modinput
    use modmain,    only : apwordmax, lmmaxapw, natmtot, nmatmax, nstfv, &
                           lmmaxvr, idxlm, nrcmt, zone, rcmt, idxas, nrcmtmax, &
                           zzero
    use mod_LDA_LU, only : llu, vmatlu, lmmaxlu
    implicit none
    integer(4), intent(in)    :: is
    integer(4), intent(in)    :: ia
    integer(4), intent(in)    :: ngk
    complex(8), intent(in)    :: apwalm(ngk,apwordmax,lmmaxapw,natmtot)
    complex(8), intent(in)    :: evecfv(nmatmax)
    complex(8), intent(inout) :: vnn
    ! local
    integer(4) :: l, nm, lm, ias
    complex(8), allocatable :: wfmt1(:,:), wfmt2(:,:)
    ! external functions
    complex(8), external :: zfmtinp

    ! compute the first-variational wavefunctions
    allocate(wfmt1(lmmaxvr,nrcmtmax))
    wfmt1 = zzero
    call wavefmt(input%groundstate%lradstep, &
                 input%groundstate%lmaxvr, is, ia, ngk, apwalm, &
                 evecfv, lmmaxvr, wfmt1)

    ias = idxas(ia,is)
    l = llu(is)
    lm = idxlm(l,-l)
    nm = 2*l+1

    allocate(wfmt2(lmmaxvr,nrcmtmax))
    wfmt2 = zzero
    call zgemm('N', 'N', nm, nrcmt(is), nm, zone, &
               vmatlu(lm,lm,1,1,ias), lmmaxlu, &
               wfmt1(lm,1), lmmaxvr, zone, wfmt2(lm,1), lmmaxvr)

    vnn = vnn + zfmtinp(.True., input%groundstate%lmaxmat, nrcmt(is), &
                        rcmt(:,is), lmmaxvr, wfmt1, wfmt2)

    deallocate(wfmt1, wfmt2)

end subroutine