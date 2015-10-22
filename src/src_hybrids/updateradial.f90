!
!BOP
! !ROUTINE: updateradial
! !INTERFACE:
!
Subroutine updateradial
! !USES:
    Use modmain
! !DESCRIPTION:
!   Update of radial functions during the Hartree-Fock hybrids run.
!
! !REVISION HISTORY:
!   Created July 2015 (UW)
!EOP
!BOC
    Implicit None

    integer :: xc_old, xcspin_old, xcgrad_old
    Integer :: is, ia, ias, ir, lm, lmmax

!----------------------------------------------------------
! create fake exact-exchange potential 
!----------------------------------------------------------
          ex_coef = (1.d0-input%groundstate%Hybrid%excoeff)
          ec_coef = 0.d0
          xc_old=xctype (1)
          xcspin_old=xcspin
          xcgrad_old=xcgrad  
          xctype (1) =  4
          xcspin=0
          xcgrad=0
          Call potxc
!----------------------------------------------------------
! add fake to effective potential  
!----------------------------------------------------------
      ! muffin-tin part
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            lmmax = lmmaxinr
            Do ir = 1, nrmt (is)
               If (ir .Gt. nrmtinr(is)) lmmax = lmmaxvr
               Do lm = 1, lmmax
                  veffmt(lm,ir,ias) = veffmt(lm,ir,ias) + vxcmt(lm,ir,ias)
               End Do
               Do lm = lmmax + 1, lmmaxvr
                  veffmt(lm,ir,ias) = 0.d0
               End Do
            End Do
         End Do
      End Do
      ! interstitial part
      veffir(:) = veffir(:) + vxcir(:)
!----------------------------------------------------------
! update radial functions  
!----------------------------------------------------------
          call genapwfr         ! generate the APW radial functions
          call genlofr(tlast)   ! generate the local-orbital radial functions
          call olprad           ! compute the overlap radial integrals
!----------------------------------------------------------
! remove fake from effective potential  
!----------------------------------------------------------
      ! muffin-tin part
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            lmmax = lmmaxinr
            Do ir = 1, nrmt (is)
               If (ir .Gt. nrmtinr(is)) lmmax = lmmaxvr
               Do lm = 1, lmmax
                  veffmt(lm,ir,ias) = veffmt(lm,ir,ias) - vxcmt(lm,ir,ias)
               End Do
               Do lm = lmmax + 1, lmmaxvr
                  veffmt(lm,ir,ias) = 0.d0
               End Do
            End Do
         End Do
      End Do
      ! interstitial part
      veffir(:) = veffir(:) - vxcir(:)
!----------------------------------------------------------
! restore original XC-definitions 
!----------------------------------------------------------
          xctype (1)=xc_old
          xcspin=xcspin_old
          xcgrad=xcgrad_old  

End Subroutine
!EOC

