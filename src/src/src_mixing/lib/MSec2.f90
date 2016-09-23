      subroutine MSEC2(Y,S,SY,F,STEP,MAXMIX,MEMORY,DMIX,IFAIL,DELTA)

           implicit real*8 (a-h,o-z)
           dimension STEP(MAXMIX),F(MAXMIX)
           dimension Y(MAXMIX,MEMORY),S(MAXMIX,MEMORY)
           dimension SY(MEMORY,MEMORY)
           real*8,allocatable :: SYINV(:,:),YPROJ(:),YTG(:), YT(:,:)
           real*8,allocatable :: Undef(:), Def(:)
!--------------------------------------------------------------------
!          MultiSecant Update, Version 2
!          This form is equivalent to a Least-Squares GB Multisecant
!
!          Input Parameters
!          MAXMIX:         Length of the variable vector
!          MEMORY:         Total number of memory values to use
!                          Most recent is last
!          Y,S:            Conventional Y and S arrays
!          SY:             Matrix of S*Y values
!          F:              -Grad(MAXMIX) at the current point (residue)
!          DMIX:           Scaler for initial matrix
!
!          Output
!          STEP            Multi-Secant Step
!
!       Delta: Term to add to Diagonal for inverse
!       Parameter (Delta=1D-3)
!
!--------------------------------------------------------------------
        allocate (SYINV(MEMORY,MEMORY),YPROJ(MEMORY),YTG(MEMORY))
        allocate (YT(MEMORY,MEMORY) )
!       Invert YY
        ninverts=0
!       Note: SY is negative here
        YT=SY
        DO N=1,MEMORY
                YT(N,N)=SY(N,N)-Delta
        ENDDO
!
        ADD=-1D-4
1000    continue
        sYINV=YT
        call minverse(sYINV,MEMORY,MEMORY,IFAIL)
!       If the matrix is singular, add a diagonal component to it
        if(ifail .lt. 0)then
                  ninverts=ninverts+1
                  do J=1,MEMORY
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
        DO N=1,MEMORY
!       Dot product Y_T . F
           YTG(N)=dot_product(S(1:MAXMIX,N),F(1:MAXMIX))
        ENDDO
!       Multiply by SYINV
        DO N=1,MEMORY
           YPROJ(N)=dot_product(YTG(1:MEMORY),SYINV(N,1:MEMORY))
        ENDDO
!       Create the step
        DO J=1,MAXMIX
!          New Position
           ST = 0
           SP = F(J)
           DO N=1,MEMORY
              ST = ST - S(J,N)*YPROJ(N)
              SP = SP - Y(J,N)*YPROJ(N)
           ENDDO
           STEP(J)  =ST+DMIX*SP
        ENDDO

        IFAIL=0
100     deallocate (SYINV, YTG, YPROJ, YT )
        RETURN
        END

