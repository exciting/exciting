subroutine setup_YY(n,S,Y,YY)
use modmixermsec,only:noldstepsmax,PWHIST,FHIST,CLMHIST,YHIST
implicit none
integer,intent(in)::n
real(8),intent(in)::S(n,noldstepsmax),Y(n,noldstepsmax)
real(8),intent(out)::YY(noldstepsmax,noldstepsmax)

real(8),parameter:: rtrap =0.1
#ifdef none
!--------------------------------------------------------------------
!       Generate the MEMORY x MEMORY Matrices
!       Also generate scaling information
!--------------------------------------------------------------------
!
!       Note: should use dgemms here -- for later
!
        DO J=1,MEMORY
          DO K=1,J
                SS(J,K) =dot_product(S (1:MAXMIX,J),S (1:MAXMIX,K))
                YY(J,K) =dot_product(Y (1:MAXMIX,J),Y (1:MAXMIX,K))
          enddo

        enddo
!       Do transpose part as well
        DO J=1,MEMORY
                DO K=1,J-1
                        SS(K,J)=SS(J,K)
                        YY(K,J)=YY(J,K)
                enddo
        enddo
         call LimitDMIX(Y,S,YY,FCURRENT,FHIST,YHIST,PWHIST,CLMHIST,MAXMIX,MEMALL,MEMORY, &
                             ascl1, qmx_input, dmix_last, qmx, dmixout, &
                             IDSCALE,Nplane, PM1, rtrap, Dbase, &
                             ISCL, RedPred, RedOld, DIAG, MUSE,MSECINFO)
#endif
end subroutine
