!
!
!
! Copyright (C) 2021 J. Cimurs, A. Gulans
!
!
Subroutine calcACE (ikp, vnlvv, vxpsiir,vxpsimt)
      Use modmain
      Use modinput
      Use modgw, only : kqset,Gkqset, kset, nomax, numin, ikvbm, ikcbm, ikvcm, Gset
      Use modmpi, only : rank
      Implicit None
! arguments
      Integer, Intent (In) :: ikp
      Complex (8), Intent (In) :: vnlvv (nstsv, nstsv)
      Complex (8), Intent (In) :: vxpsiir (ngkmax, nstsv)
      Complex (8), Intent (In) :: vxpsimt (lmmaxvr, nrcmtmax, natmtot, nstsv)

! local variables
      Integer :: ngknr, ik, jk, ist1, ist2, ist3, is4
      Integer :: is, ia, ias, ic, m1, m2, lmax, lm, ilm, irc ! ilm and lm potentially redundant (one of them)
      Integer :: nrc, iq, ig, iv (3), igq0, igk
      Integer :: ilo, loindex
      Integer :: info

      Real (8) :: v (3), cfq, ta,tb
      Complex (8) zrho01, zrho02, zt1, zt2
      Integer :: nr, l, m, io1, lm2, ir, if3
! automatic arrays
      Real (8) :: zn (nspecies)
      Complex (8) sfacgq0 (natmtot)
! allocatable arrays
      Integer, Allocatable :: igkignr (:)
      Real (8), Allocatable :: vgklnr (:, :)
      Real (8), Allocatable :: vgkcnr (:, :)
      Real (8), Allocatable :: gkcnr (:)
      Real (8), Allocatable :: tpgkcnr (:, :)
      Real (8), Allocatable :: vgqc (:, :)
      Real (8), Allocatable :: tpgqc (:, :)
      Real (8), Allocatable :: gqc (:)
      Real (8), Allocatable :: jlgqr (:, :, :)
      Real (8), Allocatable :: jlgq0r (:, :, :)
      Real (8), Allocatable :: evalsvp (:)
      Real (8), Allocatable :: evalsvnr (:)
      Real (8), Allocatable :: evalfv (:,:)
      Complex (8), Allocatable :: apwalm (:, :, :, :)
      Complex (8), Allocatable :: evecfv (:, :)
      Complex (8), Allocatable :: evecsv (:, :)
      Complex (8), Allocatable :: sfacgknr (:, :)
      Complex (8), Allocatable :: ylmgq (:, :)
      Complex (8), Allocatable :: sfacgq (:, :)
      Complex (8), Allocatable :: wfmt1 (:, :, :, :, :)
      Complex (8), Allocatable :: wfmt2 (:, :, :, :, :)
      Complex (8), Allocatable :: wfir1 (:, :, :)
      Complex (8), Allocatable :: wfir2 (:, :, :)
      Complex (8), Allocatable :: wfcr1 (:, :, :)
      Complex (8), Allocatable :: wfcr2 (:, :, :)
      Complex (8), Allocatable :: zrhomt (:, :, :)
      Complex (8), Allocatable :: zrhoir (:)
      Complex (8), Allocatable :: zvcltp (:, :)
      Complex (8), Allocatable :: zfmt (:, :)
      Complex (8), Allocatable :: matrixl(:, :)
      Complex (8), Allocatable :: matrixm(:, :)
      Complex (8), Allocatable :: matrixm1(:, :)
      Complex (8), Allocatable :: matrixm2(:, :)
      Complex (8), Allocatable :: hfxiir(:, :)
      Complex (8), Allocatable :: hfximt(:, :, :, :)
      Complex (8), Allocatable :: zfft(:)
      Complex (8) :: fr
      Real (8), Allocatable :: cf(:,:), frre (:), frim(:), gr(:)
      Real (8) :: xiintegralre, xiintegralim
      Complex (8), Allocatable :: xiintegral (:, :)
      Complex (8), Allocatable :: apwi(:, :)
      Complex (8), Allocatable :: matrixPC(:, :)
      type (WFType) :: wf1,wf2,prod,pot
! external functions
      Complex (8) zfinp, zfmtinp
      External zfinp, zfmtinp
! for decomposition
      integer :: lwork, k, len
      complex(8), allocatable :: work(:)
      complex(8), allocatable ::  diag(:,:), matrixl1(:,:)
      real(8), allocatable :: rwork(:), eigvals(:)




! -- Adaptively Compressed Exchange Operator starts --
      Allocate (matrixl(nstsv,nstsv))
      Allocate (matrixm(nstsv,nstsv))
      Allocate (matrixm1(nstsv,nstsv))
      Allocate (matrixm2(nstsv,nstsv))
      Allocate (hfxiir(ngkmax,nstsv))
      Allocate (hfximt(lmmaxvr, nrcmtmax, natmtot, nstsv))

      if (.not.(allocated(pace))) then
        Allocate (pace(nstsv, nmatmax, nkpt))
        pace=zzero
      endif
      Allocate (cf (3, nrmtmax), frre(nrmtmax), frim(nrmtmax), gr(nrmtmax))
      Allocate (xiintegral(nstsv,haaijSize))
      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate (apwi(haaijSize,ngkmax))
      Allocate (matrixPC(nstsv,nstfv))

matrixl = -vnlvv ! create copy

! Remove upper triangular part
Do ist1 = 2, nstsv
    matrixl(ist1-1,ist1:)=0
End Do


        Allocate (diag(nstsv,nstsv))
        Allocate (work(nstsv))
        Allocate (rwork(3*nstsv-2))
        Allocate (eigvals(nstsv))


        call zheev('V', 'L', nstsv, matrixl,nstsv, eigvals, work, -1, rwork, info)
        k = int(work(1))

        deallocate(work)
        allocate(work(k))

        call zheev('V', 'L', nstsv, matrixl, nstsv, eigvals, work, k, rwork, info)
        !write(*,*)info, "eigval info"


        !--sqrt matrix
        diag = 0.d0
        do k =1, nstsv
                if (eigvals(k).gt.(1e-8)) then
                        diag(k,k) = cmplx(1/sqrt(eigvals(k)),0.d0)
                else 
!        write(*,*)"small eigval", eigvals(k)
                        diag(k,k) = cmplx(0.d0, 0.d0)
                end if
        end do


        matrixl = conjg(transpose(matrixl))
        allocate(matrixl1(nstsv, nstsv))
        matrixl1 = 0.d0
        call zgemm('N', 'N', nstsv, nstsv, nstsv,    dcmplx(1.0D0,0.0), diag, nstsv, matrixl, nstsv, dcmplx(0.0D0,0.0), matrixl1, nstsv)   !construct B = sqrt(D^-1)*L^-1
        matrixl = matrixl1
        
        deallocate(work, rwork, eigvals, diag, matrixl1)



        Call zgemm('N','C',ngkmax,nstsv,nstsv,(1.0D0,0.0),vxpsiir,ngkmax,matrixl,nstsv,(0.0D0,0.0),hfxiir,ngkmax) ! compute IR part of ACEO   hfxiir=zvclir*matrixl+


Do ilm = 1, lmmaxvr
    Do irc = 1, nrcmtmax
        Call zgemm('N','C', natmtot, nstsv, nstsv, (1.0D0,0.0), vxpsimt(ilm,irc,:,:), natmtot, matrixl, nstsv, (0.0D0,0.0), hfximt(ilm,irc,:,:), natmtot) ! compute MT part of ACEO
    End Do
End Do

! -- calculating IR part of <xi|phi>
pace(:,:,ikp)=0.d0
Do ist1 = 1, nstsv
  pace(ist1,1:ngk(1,ikp), ikp) = conjg(hfxiir(1:ngk(1,ikp),ist1))
End Do

! -- calculating MT and LO part of <xi|phi>
! apwi=0
Call match (ngk(1, ikp), gkc(:, 1, ikp), tpgkc(:, :, 1, ikp), sfacgk(:, :, 1, ikp), apwalm)
loindex=0
Do is = 1, nspecies
    nr = nrmt (is)
    Do ia = 1, natoms (is)
        ias = idxas (ia, is)
        if3=0
        Do l = 0, input%groundstate%lmaxmat

                Do m = -l, l
            Do io1 = 1, apword (l, is)                
                    lm2 = idxlm (l, m)
                    ! m2 = mfromlm(lm2)
                    ! l2 = lfromlm(lm2)
                    if3=if3+1
                    apwi(if3,1:ngk(1,ikp))=apwalm(1:ngk(1,ikp), io1, lm2, ias)
                    Do ist2 = 1, nstsv
                        Do ir = 1, nr
                            fr=apwfr(ir,1,io1,l,ias)*conjg(hfximt(lm2,ir,ias,ist2)) *spr(ir, is) ** 2 ! r2(ir)=spr(ir, is) ** 2
                            frre (ir)=dble(fr)
                            frim (ir)=aimag(fr)
                        End Do
                        Call fderiv (-1, nr, spr(:, is), frre, gr, cf)
                        xiintegralre=gr (nr) ! real part
                        Call fderiv (-1, nr, spr(:, is), frim, gr, cf)
                        xiintegralim=gr (nr) ! imag part
                        xiintegral (ist2,if3)=dcmplx(xiintegralre,xiintegralim) ! nāk klāt is un ia
                    End Do
                End Do
            End Do
        End Do
! calculate LO part
        Do ilo= 1, nlorb(is)
            l = lorbl (ilo, is)
            Do m = -l, l
                    lm2 = idxlm (l, m)
                    loindex=loindex+1
                    Do ist2 = 1, nstsv
                        Do ir = 1, nr
                            fr=lofr(ir,1,ilo,ias)*conjg(hfximt(lm2,ir,ias,ist2)) *spr(ir, is) ** 2 ! r2(ir)=spr(ir, is) ** 2
                            frre (ir)=dble(fr)
                            frim (ir)=aimag(fr)
                        End Do
                        Call fderiv (-1, nr, spr(:, is), frre, gr, cf)
                        xiintegralre=gr (nr) ! real part
                        Call fderiv (-1, nr, spr(:, is), frim, gr, cf)
                        xiintegralim=gr (nr) ! imag part
!                        if(pace (ist2,loindex+ngk(1,ikp),ikp) .ne. 0.0) write(*,*) 'paceLO not zero'
                        pace (ist2,loindex+ngk(1,ikp),ikp)= pace (ist2,loindex+ngk(1,ikp),ikp)+dcmplx(xiintegralre,xiintegralim) ! paceLO
                        !write(*,*) gr(1), gr(nr), gr(150), gr(151)
                    End Do
            End Do
        End Do ! LO
        Call zgemm('N', 'N', nstsv, ngk(1,ikp), if3, dcmplx(1.0D0,0.0), xiintegral, nstsv, apwi, haaijSize, dcmplx(1.0D0,0.0), pace(:,:,ikp), nstsv) ! pace= paceMT+pace(IR+LO) = xiintegral*apwi+ pace
    End Do
End Do

!! -- Adaptively Compressed Exchange Operator test
if (.false.) then
if (rank==0) then
    Allocate (evecfv(nmatmax, nstfv))
    Call getevecfv (vkl(:, ikp), vgkl(:, :, :, ikp), evecfv)
    call WFInit(wf1)
    call genWF(ikp,wf1)
    call genWFinMT(wf1)    
! test pace
    Call zgemm('N', 'N', nstsv, nstfv, nmat(1,ikp), dcmplx(1.0D0,0.0), pace(:,:,ikp), nstsv, evecfv, nmatmax, dcmplx(0.0D0,0.0), matrixPC, nstsv) ! matrixPC=pace*evecfv
    Call zgemm('C', 'N', nstfv, nstfv, nstsv, dcmplx(-1.0D0,0.0), matrixPC, nstfv, matrixPC, nstsv, dcmplx(0.0D0,0.0), matrixm2, nstfv) ! matrixM2=-matrixPC^+ * matrixPC

! test hfxi
    Do ist1 = 1, nstsv
        Do ist2 = 1, nstsv
            matrixm1(ist1,ist2) = zfinp(.True., wf1%mtrlm(:,:,:,ist1), hfximt(:,:,:,ist2), wf1%ir(:,ist1), hfxiir(:,ist2)) ! matrixM1=<phi|hfxi>
        End Do
    End Do
    Call zgemm('N','C', nstsv, nstsv, nstsv, dcmplx(-1.0D0,0.0), matrixm1, nstsv, matrixm1, nstsv, dcmplx(0.0D0,0.0), matrixm, nstsv) ! matrixM=-matrixM1 * matrixM1^+

! matrix print compare
    write(*,*) 'ikp=',ikp
    write(*,*) 'matrixm2 real (pace)'
      do ist1 = 1, 12
        write(*,'(12F13.9)') dble(matrixm2(ist1,1:12))
      end do

    write(*,*) 'matrixm1 real (hfxi)'
      do ist1 = 1, 12
        write(*,'(12F13.9)') dble(matrixm(ist1,1:12))
      end do

    write(*,*) 'matrixm real (vnlvv)'
      do ist1 = 1, 12
        write(*,'(12F13.9)') dble(vnlvv(ist1,1:12))
      end do

    write(*,*) 'matrixm2 imag (pace)'
      do ist1 = 1, 12
        write(*,'(12F13.9)') dimag(matrixm2(ist1,1:12))
      end do

    write(*,*) 'matrixm2 imag (hfxi)'
      do ist1 = 1, 12
        write(*,'(12F13.9)') dimag(matrixm(ist1,1:12))
      end do

    write(*,*) 'matrixm imag (vnlvv)'
      do ist1 = 1, 12
        write(*,'(12F13.9)') dimag(vnlvv(ist1,1:12))
      end do
    write(*,*) '--------------'
     deallocate(evecfv)
     call WFRelease(wf1)
endif
endif

  Deallocate (matrixl, matrixm, matrixm1, matrixm2)
  Deallocate (hfxiir,hfximt)
  Deallocate(cf,frre,frim,gr)
  Deallocate(xiintegral)
  Deallocate (apwi,apwalm)
  Deallocate (matrixPC)


!-- Adaptively Compressed Exchange Operator ends --

      Return
End Subroutine
!EOC
