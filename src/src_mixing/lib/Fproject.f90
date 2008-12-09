      real(8) function Fprojmem(Y,S,F,YY,MAXMIX,MEMORY,PSTEP, &
      DELTA,NPLANE,MUSE)
!     Projects the current residue onto prior gradient information
      implicit real*8 (a-h,o-z)
      dimension Y(MAXMIX,MEMORY), YY(MEMORY,MEMORY), S(MAXMIX,MEMORY)
      dimension F(MAXMIX)
      real*8,allocatable :: YYINV(:,:),YPROJ(:),YTG(:), YT(:,:)
      real*8,allocatable :: ETA(:), XI(:) !, TTT(:,:)
!--------------------------------------------------------------------
!       Estimate the NFL Projection
!
!       Input Parameters
!       MAXMIX:         Length of the variable vector
!       MEMORY:         Total number of memory values to use
!                          Most recent is last
!       Y:              Conventional Y arrays
!       YY:             Matrix of Y*Y values
!       F:              -Grad(MAXMIX) at the current point (residue)
!
!       Output
!       FProject        Ratio of Undefined/Defined
!
!       Delta: Term to add to Diagonal for inverse
!       Parameter (Delta=5D-4)
!
!--------------------------------------------------------------------
      allocate (YYINV(MUSE,MUSE),YPROJ(MUSE),YTG(MUSE))
      allocate (YT(MUSE,MUSE) )
      allocate (ETA(MAXMIX), XI(MAXMIX) )
        fprojmem=-1.D0

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
!                 Bomb
                  goto 100
        endif
!
!       Mini Hack
!        allocate (TTT(MUSE,MUSE))
!        YT(1:MUSE,1:MUSE)=YY(1+ISKIP:MEMORY,1+ISKIP:MEMORY)
!        DO N=1,MUSE
!         DO L=1,MUSE
!           T=0
!           DO M=1,MUSE
!                T=T+YYINV(N,M)*YT(M,L)
!           ENDDO
!           TTT(N,L)=T
!         ENDDO
!        ENDDO
!        T=0
!        DO N=1,MUSE
!12              format(10F12.5)
!                write(21,12)TTT(1:MUSE,N)
!                DO M=1,MUSE
!                        T=T+TTT(N,M)*TTT(M,N)
!                ENDDO
!        ENDDO
!        T=T/MUSE
!        write(21,*)'Debug ',T,1+Delta
!
        DO N=1,MUSE
!               Dot product Y_T . F
                YTG(N)=dot_product(Y(1:MAXMIX,N+ISKIP),F(1:MAXMIX))
        ENDDO
!       Multiply by YYINV
        DO N=1,MUSE
           YPROJ(N)=dot_product(YTG(1:MUSE),YYINV(1:MUSE,N))
        ENDDO
!       Create the undefined part
        Residue=0.0
        Defined=0.0
        DO J=1,MAXMIX
           SP = F(J)
           Defined = Defined+SP*SP
           DO N=1,MUSE
              SP = SP - Y(J,N+ISKIP)*YPROJ(N)
           ENDDO
           Residue = Residue + SP*SP
        ENDDO

        if(Residue .gt. 1D-30)then
                FProjmem=sqrt(Residue/Defined)
        else
                FProjmem=0.0D0
        endif
!       Here Eta is the projected known gradient, XI is the projected known step
        ETA = 0
        XI  = 0
        DO N=1,MUSE
        ETA(1:MAXMIX) = ETA(1:MAXMIX)+YPROJ(N)*Y(1:MAXMIX,N+ISKIP)
        XI (1:MAXMIX) = XI (1:MAXMIX)-YPROJ(N)*S(1:MAXMIX,N+ISKIP)
        ENDDO
        EE = dot_product(ETA,ETA)
        EX = dot_product(ETA,XI )
        PSTEP = abs(EX/EE)

100     deallocate(YT, YTG, YPROJ, YYINV,ETA,XI) !,TTT)
        return
        end

           real*8 Function BBProject(Y,S,YY,F,MAXMIX,MEMORY)
           implicit real*8 (a-h,o-z)
           dimension F(MAXMIX),YY(MEMORY,MEMORY)
           dimension Y(MAXMIX,MEMORY),S(MAXMIX,MEMORY)
           real*8,allocatable :: WTMP(:),ZTMP(:)
!--------------------------------------------------------------------
!          Find the BB residue to known ratio
!
!          Input Parameters
!          MAXMIX:         Length of the variable vector
!          MEMORY:         Total number of memory values to use
!                          Most recent is last
!          Y,S:            Conventional Y and S arrays
!          YY:             Matrix of Y*Y values
!          F:              -Grad(MAXMIX) at the current point (residue)
!
!          Output
!          Value           |Unknown|/|Input|
!--------------------------------------------------------------------
!
           allocate(WTMP(MAXMIX),ZTMP(MAXMIX))
           TEMPZW= -DOT_PRODUCT(Y(1:MAXMIX,MEMORY),F(1:MAXMIX))
           TEMPZW=  TEMPZW/YY(MEMORY,MEMORY)

           Wtmp(1:MAXMIX)= -F(1:MAXMIX)
           Ztmp(1:MAXMIX)=  S(1:MAXMIX,MEMORY)*TEMPZW
           DO J=0,MEMORY-2
              TEMPW=TEMPZW
              TEMPZW=0.0D0
              DO K=1,MAXMIX
                 Wtmp(K) = Wtmp(K)- Y(K,MEMORY-J)*TEMPW
                 TEMPZW =  TEMPZW+Y(K,MEMORY-J-1)*Wtmp(K)
              ENDDO
              TEMPZW = TEMPZW/YY(MEMORY-J-1,MEMORY-J-1)
              Ztmp(1:MAXMIX)=Ztmp(1:MAXMIX)+S(1:MAXMIX,MEMORY-J-1)*TEMPZW
           ENDDO
!
!           Finish up
           Wtmp(1:MAXMIX)=Wtmp(1:MAXMIX) - Y(1:MAXMIX,1)*TEMPZW
!
           Defined=dot_product(F(1:maxmix),F(1:maxmix))
!          Residue
           Residue=dot_product(Wtmp(1:maxmix),Wtmp(1:maxmix))
           Defined=sqrt(Defined)
           if(Residue .gt. 1D-30)then
                BBProject=sqrt(Residue)/Defined
           else
                BBProject=0.0D0
           endif
           deallocate(Wtmp, Ztmp)

           return
           end

      Real*8 Function GBMProject(Y,S,SY,F,MAXMIX,MEMORY)
           implicit real*8 (a-h,o-z)
           dimension F(MAXMIX)
           dimension Y(MAXMIX,MEMORY),S(MAXMIX,MEMORY)
           dimension SY(MEMORY,MEMORY)
           real*8,allocatable :: SYINV(:,:),YPROJ(:),YTG(:), YT(:,:)
           real*8,allocatable :: Undef(:), Def(:)
!--------------------------------------------------------------------
!          Find the Good-Broyden (Closest Point) Residue
!
!          Input Parameters
!          MAXMIX:         Length of the variable vector
!          MEMORY:         Total number of memory values to use
!                          Most recent is last
!          Y,S:            Conventional Y and S arrays
!          SY:             Matrix of S*Y values
!          F:              -Grad(MAXMIX) at the current point (residue)
!
!          Output
!          Ration of Unknown/Known
!
!       Delta: Term to add to Diagonal for inverse
        Parameter (Delta=1D-3)
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
                  GMBProject=-999
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
!       Create the undefined part
        Residue=0.0
        Defined=0.0
        DO J=1,MAXMIX
           SP = F(J)
           Defined = Defined+SP*SP
           DO N=1,MEMORY
              SP = SP - Y(J,N)*YPROJ(N)
           ENDDO
           Residue = Residue + SP*SP
        ENDDO
        Defined=sqrt(Defined)
        if(Residue .gt. 1D-30)then
                GBMProject=sqrt(Residue)/Defined
        else
                GBMProject=0.0D0
        endif
100     deallocate (SYINV, YTG, YPROJ, YT )
        RETURN
        END

        real*8 function probegb(S,Y,SY,SS,F,MAXMIX,MEMORY,rtrap)
        implicit real*8 (a-h,o-z)
        dimension F(MAXMIX),SS(MEMORY,MEMORY)
        dimension Y(MAXMIX,MEMORY),S(MAXMIX,MEMORY),SY(MEMORY,MEMORY)
        real*8,allocatable :: SB(:)
!
!       if(MEMORY .lt. 2)return
        allocate (SB(MAXMIX))
!       Generate the GB step
        DO DMIXM=0.005,2.0,0.005
        call  GBLDM(Y,S,SY,SS,F,SB,MAXMIX,MEMORY,DMIXM,IFAIL)
!       Find the corresponding back projection, with everything flipped
        T1 = BBProject(S,Y,SS,SB,MAXMIX,MEMORY)
        WeWant = dmixm*T1/rtrap
!        write(21,*)'Test ',DMIXM,T1,(dmixm*T1/rtrap)
        if(WeWant .ge. 1.0)then
                ProbeGB=rtrap/dmixm
                deallocate(SB)
                return
                endif
        enddo
        deallocate (SB)
        call exit(0)
        end

