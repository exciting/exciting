      subroutine minverse ( A, NP , NDIM, IFAIL )
      !-------------------------------------------------------------------------
      !
      !	      Taken from "Numeric recipes".  The original program was
      !       GAUSSJ which solves linear equations by the Gauss_Jordon
      !       elimination method.  Only the parts required to invert
      !	      matrices have been retained.
      !
      !	      J.P. Griffith  6/88
      !
      !-------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      PARAMETER (NMAX=500)
      dimension A(NDIM,NDIM)
      integer   IPIV(NMAX), INDXR(NMAX), INDXC(NMAX)

      n = np
      IFAIL = 0
!     Huh, why is this needed ?
!      ICOL=1
!      irow=1
      ICOL=0
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.d0
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                IFAIL=-1
                RETURN
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        if(icol.eq.0)then
           IFAIL=-1
!           write(6,*)'Singular ?'
!           write(6,*)(A(J,K),K=1,N)
           return
        endif
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.D0) THEN
           IFAIL = -2
           RETURN
        ENDIF
        PIVINV=1.d0/A(ICOL,ICOL)
        A(ICOL,ICOL)=1.d0
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.d0
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END
