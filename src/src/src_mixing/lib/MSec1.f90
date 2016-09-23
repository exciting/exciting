     	subroutine MSEC1(Y,S,YY,F,STEP,MAXMIX,MEMORY,DMIX,IFAIL,DELTA,MUSE)
           implicit real*8 (a-h,o-z)
           dimension STEP(MAXMIX),F(MAXMIX)
           dimension Y(MAXMIX,MEMORY),S(MAXMIX,MEMORY)
           dimension YY(MEMORY,MEMORY)
           real*8,allocatable :: YYINV(:,:),YPROJ(:),YTG(:), YT(:,:)
!--------------------------------------------------------------------
!          No Free Lunch MultiSecant Update
!
!          Input Parameters
!          MAXMIX:         Length of the variable vector
!          MEMORY:         Total number of memory values to use
!                          Most recent is last
!          Y,S:            Conventional Y and S arrays
!          YY:             Matrix of Y*Y values
!          F:              -Grad(MAXMIX) at the current point (residue)
!          HZERO:          Diagonal initial matrix
!          DMIX:           Scaler for initial matrix
!
!          Output
!          STEP            Multi-Secant Step
!
!       Delta: Term to add to Diagonal for inverse
!       Parameter (Delta=1D-3)
!
!--------------------------------------------------------------------
        allocate (YYINV(MUSE,MUSE),YPROJ(MUSE),YTG(MUSE))
        allocate (YT(MUSE,MUSE) )
        ifail=1
!       Invert YY
        ninverts=0
        ISKIP=MEMORY-MUSE
        YT(1:MUSE,1:MUSE)=YY(1+ISKIP:MEMORY,1+ISKIP:MEMORY)
        DO N=1,MUSE
                YT(N,N)=YT(N,N)+Delta
        ENDDO
!
        ADD=1D-4
1000    continue
        YYINV=YT
        call minverse(YYINV,MUSE,MUSE,IFAIL)
!       If the matrix is singular, add a diagonal component to it
        if(ifail .lt. 0)then
                  ninverts=ninverts+1
                  do J=1,MUSE
                     YT(J,J)=YT(J,J)+ADD
                  enddo
                  ADD=ADD*2.D0
                  if(ninverts .lt. 10)goto 1000
!                 Bomb to Pratt
                  IFAIL=-1
                  goto 100
        endif
!
        ntries=0
2000    continue
        DO N=1,MUSE
!       Dot product Y_T . F
           YTG(N)=dot_product(Y(1:MAXMIX,N+ISKIP),F(1:MAXMIX))
        ENDDO
!       Multiply by YYINV
        DO N=1,MUSE
           YPROJ(N)=dot_product(YTG(1:MUSE),YYINV(1:MUSE,N))
        ENDDO
!       Create the step
        DO J=1,MAXMIX
!          New Position
           ST = 0
           SP = F(J)
           DO N=1,MUSE
              ST = ST - S(J,N+ISKIP)*YPROJ(N)
              SP = SP - Y(J,N+ISKIP)*YPROJ(N)
           ENDDO
           STEP(J)  =ST+DMIX*SP

        ENDDO

!   write(*,*)"DMIX IN MSEC1:",DMIX

        IFAIL=0
100     deallocate (YYINV, YTG, YPROJ, YT )
        RETURN
        END

