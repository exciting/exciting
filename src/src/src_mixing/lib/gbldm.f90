       subroutine GBLDM(Y,S,SY,SS,F,STEP,MAXMIX,MEMORY,DMIX,IFAIL)
           implicit real*8 (a-h,o-z)
           dimension STEP(MAXMIX),F(MAXMIX)
           dimension Y(MAXMIX,MEMORY),S(MAXMIX,MEMORY)
           dimension SS(MEMORY,MEMORY),SY(MEMORY,MEMORY)
           real*8,allocatable :: YYINV(:,:),YPROJ(:),YTG(:), YT(:,:)
           real*8,allocatable :: Undef(:), Def(:)
!--------------------------------------------------------------------
!          Good Broyden according to LDM
!
!          Input Parameters
!          MAXMIX:         Length of the variable vector
!          MEMORY:         Total number of memory values to use
!                          Most recent is last
!          Y,S:            Conventional Y and S arrays
!          YY:             Matrix of Y*Y values
!          SY:             Matrix of S*Y values
!          F:              -Grad(MAXMIX) at the current point (residue)
!          DMIX:           Scaler for initial matrix
!
!          Output
!          STEP            Multi-Secant Step
!
!
!--------------------------------------------------------------------
        allocate (YYINV(MEMORY,MEMORY),YPROJ(MEMORY),YTG(MEMORY))
        allocate (YT(MEMORY,MEMORY) )
!       Invert YY
        ninverts=0
!
        YT=SY
!       Difference between BB and Multi-Secant, from Byrd Nocedal Schnabel 6.2
        DO I=2,MEMORY
                DO J=1,I-1
                        YT(I,J)=YT(I,J)-SS(I,J)/dmix
                ENDDO
        ENDDO

!
!       Note, Add is negative since matrix is negative
        ADD=-1D-4
1000    continue
        YYINV=YT
        call minverse(YYINV,MEMORY,MEMORY,IFAIL)
!
!       If the matrix is singular, add a diagonal component to it
        if(ifail .lt. 0)then
                  Write(*,*)':DIRE Singular, adding ',add
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
!       Multiply by YYINV
        DO N=1,MEMORY
           YPROJ(N)=dot_product(YTG(1:MEMORY),YYINV(N,1:MEMORY))
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
100     deallocate (YYINV, YTG, YPROJ )
        RETURN
        END

