Module modmixer
  Use precision, only: dp
  real(8), allocatable :: history_ir(:, :), residual_history_ir(:, :)
  real(8), allocatable :: history_mt(:, :, :, :), residual_history_mt(:, :, :, :)
  integer :: history_size
  integer :: npsd
  logical :: densitymixing
  integer :: lmaxmix, lmmaxmix
  real(8), allocatable :: beta_ir(:), beta_mt(:, :, :)
  real(8), allocatable :: Amix(:, :)

Contains

  !> Select and execute the active mixer on rank 0, then broadcast the mixed fields.
  Subroutine pickmixer(fn_ir, fn_mt, sclstep)
    Use modinput, only: input
    Use modmain, only: natmtot, ngrtot, nrmtmax
    use exciting_mpi, only: xmpi_bcast
    use modmpi, only: mpiglobal, rank
    implicit none
    !> Current self-consistent-field iteration index.
    integer, intent(in) :: sclstep
    !> Interstitial quantity to be mixed and broadcast.
    real(8), intent(inout) :: fn_ir(ngrtot)
    !> Muffin-tin quantity to be mixed and broadcast.
    real(8), intent(inout) :: fn_mt(lmmaxmix, nrmtmax, natmtot)

    if (rank .eq. 0) then
      if (input%groundstate%mixer .eq. 'lin') then
        call adaptivemixer(fn_ir, fn_mt, sclstep)
      elseif (input%groundstate%mixer .eq. 'kerker') then
        call kerkermixer(fn_ir, fn_mt, sclstep)
      elseif (input%groundstate%mixer .eq. 'simplelinear') then
        call linearmixer(fn_ir, fn_mt, sclstep)
      elseif (input%groundstate%mixer .eq. 'pulay') then
        call pulaymixer(fn_ir, fn_mt, sclstep)
      end if
    end if

    call xmpi_bcast(mpiglobal, fn_ir)
    call xmpi_bcast(mpiglobal, fn_mt)
  End Subroutine

  !> Apply a fixed linear mixing step using the most recent residual.
  Subroutine linearmixer(fn_ir, fn_mt, sclstep)
    Use modinput, only: input
    Use modmain, only: currentconvergence, natmtot, ngrtot, nrmtmax, omega
    implicit none
    !> Interstitial quantity to be updated by the mixer.
    real(8), intent(inout) :: fn_ir(ngrtot)
    !> Muffin-tin quantity to be updated by the mixer.
    real(8), intent(inout) :: fn_mt(lmmaxmix, nrmtmax, natmtot)
    !> Current self-consistent-field iteration index.
    integer, intent(in) :: sclstep
    integer :: prevrec, nxtrec
    real(8), external :: rfinp

    prevrec = mod(sclstep - 1, history_size) + 1
    nxtrec = mod(sclstep, history_size) + 1

    residual_history_ir(:, prevrec) = fn_ir(:) - history_ir(:, prevrec)
    residual_history_mt(:, :, :, prevrec) = fn_mt(:, :, :) - history_mt(:, :, :, prevrec)

    currentconvergence = sqrt(rfinp(1, residual_history_mt(:, :, :, prevrec), residual_history_mt(:, :, :, prevrec), &
   &   residual_history_ir(:, prevrec), residual_history_ir(:, prevrec)) / omega)

    fn_ir = history_ir(:, prevrec) + input%groundstate%beta0 * residual_history_ir(:, prevrec)
    fn_mt = history_mt(:, :, :, prevrec) + input%groundstate%beta0 * residual_history_mt(:, :, :, prevrec)

    history_ir(:, nxtrec) = fn_ir
    history_mt(:, :, :, nxtrec) = fn_mt

  End Subroutine

  !> Apply Kerker screening to the residual before the linear update.
  Subroutine kerkermixer(fn_ir, fn_mt, sclstep)
        Use modinput, only: input
        Use modmain, only: currentconvergence, fourpi, natmtot, ngrtot, nrmtmax, omega
        Use modmixer_screening, Only: potcorr
        implicit none
        !> Interstitial quantity to be updated by the mixer.
        real(8), intent(inout) :: fn_ir(ngrtot)
        !> Muffin-tin quantity to be updated by the mixer.
        real(8), intent(inout) :: fn_mt(lmmaxmix, nrmtmax, natmtot)
        !> Current self-consistent-field iteration index.
        integer, intent(in) :: sclstep
        real(8), allocatable :: screened_ir(:)
        real(8), allocatable :: screened_mt(:, :, :)
        integer :: prevrec, nxtrec
        real(8), external :: rfinp


        prevrec = mod(sclstep - 1, history_size) + 1
        nxtrec = mod(sclstep, history_size) + 1

        residual_history_ir(:, prevrec) = fn_ir(:) - history_ir(:, prevrec)
        residual_history_mt(:, :, :, prevrec) = fn_mt(:, :, :) - history_mt(:, :, :, prevrec)

        currentconvergence = sqrt(rfinp(1, residual_history_mt(:, :, :, prevrec), residual_history_mt(:, :, :, prevrec), &
       &   residual_history_ir(:, prevrec), residual_history_ir(:, prevrec)) / omega)

        allocate(screened_ir(ngrtot))
        allocate(screened_mt(lmmaxmix, nrmtmax, natmtot))
        screened_ir = 0d0
        screened_mt = 0d0

        call potcorr(residual_history_mt(:, :, :, prevrec), residual_history_ir(:, prevrec), input%groundstate%lambda, &
       &   screened_mt, screened_ir)

        fn_ir = history_ir(:, prevrec) + input%groundstate%beta0 * residual_history_ir(:, prevrec)
        fn_ir = fn_ir - ((input%groundstate%beta0 * (input%groundstate%lambda**2)) / fourpi) * screened_ir
        fn_mt = history_mt(:, :, :, prevrec) + input%groundstate%beta0 * residual_history_mt(:, :, :, prevrec)
        fn_mt = fn_mt - ((input%groundstate%beta0 * (input%groundstate%lambda**2)) / fourpi) * screened_mt

        deallocate(screened_ir, screened_mt)
        
        history_ir(:, nxtrec) = fn_ir
        history_mt(:, :, :, nxtrec) = fn_mt

  End Subroutine

      !> Build a Pulay/Broyden-style mixed field from the stored residual history.
      Subroutine pulaymixer(fn_ir,fn_mt,sclstep)
        Use modinput, only: input
        Use modmain, only: currentconvergence, efermi, fermidos, fourpi, natmtot, ngrtot, nrmtmax, &
     &    omega
        Use modmpi, only: terminate_if_false
        Use modmixer_screening, Only: potcorr
        implicit none
        !> Interstitial quantity to be updated by the mixer.
        real(8),intent(inout) :: fn_ir(ngrtot)
        !> Muffin-tin quantity to be updated by the mixer.
        real(8),intent(inout) :: fn_mt(lmmaxmix,nrmtmax,natmtot)
        !> Current self-consistent-field iteration index.
        integer,intent(in) :: sclstep
        integer :: prevrec,nxtrec,maxrec
        integer :: irec,jrec
        integer :: info
        real(8), allocatable :: Amixinv(:,:),work(:),evals(:)
        integer, allocatable :: ipiv(:)
        real(8) :: factor
        real(8), allocatable :: transformed_ir(:)
        real(8), allocatable :: transformed_mt(:,:,:)
        real(8),external :: rfinp


        allocate(transformed_ir(ngrtot))
        allocate(transformed_mt(lmmaxmix,nrmtmax,natmtot))

        transformed_mt=0d0
        transformed_ir=0d0

        prevrec=mod(sclstep-1,history_size)+1
        nxtrec=mod(sclstep,history_size)+1
        
        residual_history_ir(:,prevrec)=fn_ir(:)-history_ir(:,prevrec)
        residual_history_mt(:,:,:,prevrec)=fn_mt(:,:,:)-history_mt(:,:,:,prevrec)

        if (sclstep.lt.history_size) then
          maxrec=prevrec
        else  
          maxrec=history_size
        endif

        do irec=1,maxrec
          do jrec=irec,maxrec
             call residualproduct(residual_history_ir(:,irec),residual_history_mt(:,:,:,irec),residual_history_ir(:,jrec),residual_history_mt(:,:,:,jrec),Amix(irec,jrec))
          end do
        end do

       allocate(Amixinv(maxrec+1,maxrec+1),ipiv(maxrec+1),evals(maxrec+1),work(maxrec+1))

        Amixinv=1d0
        do irec=1,maxrec
          do jrec=1,maxrec
             Amixinv(irec,jrec)=Amix(irec,jrec)
             Amixinv(maxrec+1,maxrec+1)=0d0
          enddo
        enddo
   
        evals(:)=0d0
        evals(maxrec+1)=1d0

        call dsysv('U',maxrec+1,1, Amixinv, maxrec+1,ipiv, evals, maxrec+1,work,maxrec+1,info)   
        call terminate_if_false(info .eq. 0, &
     &   'Error(pulaymixer): dsysv failed while solving the Pulay mixing system')

        fn_ir=0d0
        fn_mt=0d0

        if (sclstep.gt.input%groundstate%PrePKSteps) Then
          do irec=1,maxrec
            fn_ir=fn_ir+evals(irec)*(history_ir(:,irec)+(input%groundstate%beta0*residual_history_ir(:,irec)))
            fn_mt=fn_mt+evals(irec)*(history_mt(:,:,:,irec)+(input%groundstate%beta0*residual_history_mt(:,:,:,irec)))
          enddo
        else
          factor=max(sqrt(fourpi*abs(fermidos/omega*efermi)),input%groundstate%lambda)
          do irec=1,maxrec
            fn_ir=fn_ir+evals(irec)*(residual_history_ir(:,irec))
            fn_mt=fn_mt+evals(irec)*(residual_history_mt(:,:,:,irec))
          enddo
          call potcorr(fn_mt(:,:,:), fn_ir(:), factor, transformed_mt, transformed_ir)
          fn_ir=-(input%groundstate%beta0)*((factor**2)/fourpi)*transformed_ir
          fn_mt=-(input%groundstate%beta0)*((factor**2)/fourpi)*transformed_mt

          do irec=1,maxrec
            fn_ir=fn_ir+evals(irec)*(history_ir(:,irec)+(input%groundstate%beta0*residual_history_ir(:,irec)))
            fn_mt=fn_mt+evals(irec)*(history_mt(:,:,:,irec)+(input%groundstate%beta0*residual_history_mt(:,:,:,irec)))
          enddo
        endif

        history_ir(:,nxtrec)=fn_ir
        history_mt(:,:,:,nxtrec)=fn_mt
        currentconvergence=sqrt(rfinp(1,residual_history_mt(:,:,:,prevrec),residual_history_mt(:,:,:,prevrec),residual_history_ir(:,prevrec),residual_history_ir(:,prevrec))/omega)
        deallocate(work,ipiv,evals)
        deallocate(Amixinv)
        deallocate(transformed_ir,transformed_mt)
   
     End Subroutine

     !> Adapt the local mixing factor from the sign change history of the residual.
     Subroutine adaptivemixer(fn_ir,fn_mt,sclstep)
        Use modinput, only: input
        Use modmain, only: currentconvergence, idxas, natmtot, natoms, ngrtot, nrmt, nrmtmax, &
     &    nspecies, omega
        implicit none
        !> Interstitial quantity to be updated by the mixer.
        real(8),intent(inout) :: fn_ir(ngrtot)
        !> Muffin-tin quantity to be updated by the mixer.
        real(8),intent(inout) :: fn_mt(lmmaxmix,nrmtmax,natmtot)
        !> Current self-consistent-field iteration index.
        integer,intent(in) :: sclstep
        integer :: prevrec,nxtrec,oldrec
        integer :: ia,is,ias,ir,lm
        real(8) :: t
        real(8),external :: rfinp
 
        oldrec=mod(sclstep-2,history_size)+1
        prevrec=mod(sclstep-1,history_size)+1
        nxtrec=mod(sclstep,history_size)+1

        do ir=1,ngrtot
          t=fn_ir(ir)-history_ir(ir,prevrec)
          if (t*residual_history_ir(ir,oldrec).gt.0d0) then
            beta_ir(ir)=beta_ir(ir)*input%groundstate%betainc
            if (beta_ir(ir).gt.1d0) beta_ir(ir)=1d0
          else
            beta_ir(ir)=beta_ir(ir)*input%groundstate%betadec
          endif
          residual_history_ir(ir,prevrec)=t
        enddo

        do is=1,nspecies
          do ia=1,natoms(is)
            ias=idxas(ia,is)
            do ir=1,nrmt(is)
              do lm=1,lmmaxmix
                t=fn_mt(lm,ir,ias)-history_mt(lm,ir,ias,prevrec)
                if (t*residual_history_mt(lm,ir,ias,oldrec).gt.0d0) then
                  beta_mt(lm,ir,ias)=beta_mt(lm,ir,ias)*input%groundstate%betainc
                  if (beta_mt(lm,ir,ias).gt.1d0) beta_mt(lm,ir,ias)=1d0 
                else
                  beta_mt(lm,ir,ias)=beta_mt(lm,ir,ias)*input%groundstate%betadec 
                endif
                residual_history_mt(lm,ir,ias,prevrec)=t
              enddo
            enddo
          enddo
        enddo

        currentconvergence=sqrt(rfinp(1,residual_history_mt(:,:,:,prevrec),residual_history_mt(:,:,:,prevrec),residual_history_ir(:,prevrec),residual_history_ir(:,prevrec))/omega)

        do ir=1,ngrtot
          fn_ir(ir)=history_ir(ir,prevrec)+beta_ir(ir)*residual_history_ir(ir,prevrec)
          history_ir(ir,nxtrec)=fn_ir(ir)
        enddo

        do is=1,nspecies
          do ia=1,natoms(is)
            ias=idxas(ia,is)
            do ir=1,nrmt(is)
              do lm=1,lmmaxmix
                fn_mt(lm,ir,ias)=history_mt(lm,ir,ias,prevrec)+beta_mt(lm,ir,ias)*residual_history_mt(lm,ir,ias,prevrec)
                history_mt(lm,ir,ias,nxtrec)=fn_mt(lm,ir,ias)
              enddo
            enddo
          enddo
        enddo

      End Subroutine

      !> Evaluate the residual inner product used to build the Pulay mixing matrix.
      Subroutine residualproduct(resA_ir,resA_mt,resB_ir,resB_mt,answer)
        Use modinput, only: input
        Use modmain, only: cfunir, fourpi, gc, idxas, natmtot, natoms, ngrtot, nrmt, nrmtmax, &
     &    nspecies, omega, spr
        Use modmixer_coulomb, only: potcoulmixer
        implicit none
        !> First interstitial residual vector.
        real(8),intent(in) :: resA_ir(ngrtot)
        !> First muffin-tin residual field.
        real(8),intent(in) :: resA_mt(lmmaxmix,nrmtmax,natmtot)
        !> Second interstitial residual vector.
        real(8),intent(in) :: resB_ir(ngrtot)
        !> Second muffin-tin residual field.
        real(8),intent(in) :: resB_mt(lmmaxmix,nrmtmax,natmtot)
        !> Residual inner product value.
        real(8),intent(out) :: answer
        real(8) :: fr(nrmtmax),cf(3,nrmtmax),gr(nrmtmax)
        real(8) :: resmtrc_mt(lmmaxmix,nrmtmax,natmtot),resmtrc_ir(ngrtot)
        real(8) :: summa
        integer :: ia,is,ias,ir,lm

        summa=0d0
        if (input%groundstate%mixerswitch.eq.1) then
          if (.false.) then
            summa = sum(resA_ir(:) * resB_ir(:) * cfunir(:))
            summa = omega * summa / real(ngrtot, dp)
          endif

          do is=1,nspecies
            do ia=1,natoms(is)
              ias=idxas(ia,is)
              do lm=1,lmmaxmix
                do ir=1,nrmt(is)
                  fr(ir)=spr(ir, is)**2*resA_mt(lm,ir,ias)*resB_mt(lm,ir,ias)
                enddo
                Call fderiv (-1, nrmt(is), spr(:,is), fr, gr, cf)
                summa=summa+gr(nrmt(is))
              enddo
            enddo
          enddo
        else
          call potcoulmixer(resB_mt,resB_ir,resmtrc_mt,resmtrc_ir)
          if (.true.) then
            summa = sum(((resA_ir(:) * resB_ir(:)) + &
           &  (resA_ir(:) * (((20d0 * gc(2))**2) / fourpi) * resmtrc_ir(:))) * cfunir(:))
            summa = omega * summa / real(ngrtot, dp)
          endif

          do is=1,nspecies
            do ia=1,natoms(is)
               ias=idxas(ia,is)
              do lm=1,lmmaxmix
                do ir=1,nrmt(is)
                  fr(ir)=((resA_mt(lm,ir,ias)*resB_mt(lm,ir,ias))+(resA_mt(lm,ir,ias)*(((20d0*gc(2))**2)/fourpi)*resmtrc_mt(lm,ir,ias)))*(spr(ir, is)**2)
                enddo
                Call fderiv (-1, nrmt(is), spr(:,is), fr, gr, cf)
                summa=summa+gr(nrmt(is))
              enddo
            enddo
          enddo
        end if
        answer=summa 

      End Subroutine  

End Module
