        SUBROUTINE qmix7(clmnew,clmold,roknew,rokold,jspin,nq3,nat,&
                 lmmax,jri,qmx, ndm,jatom1,ll,alx,aly,alz,dmat,alx_old,&
                 aly_old,alz_old,dmat_old,mult,inversion,scl1,scl2,&
                 nbroyd,nmax,scale,LM,ndif, KVECVL,mode,VI,MIX,&
                 qmx_input,qtot,zz,fname, usererror)
        IMPLICIT REAL*8 (A-H,O-Z)
        INCLUDE 'param.inc'
!
        INTEGER, INTENT(IN)       :: ndm,jatom1(nat,4,3),ll(nat,4,3)
        REAL*8, INTENT(IN)        :: alx_old(nat,4,3),aly_old(nat,4,3),alz_old(nat,4,3)
        COMPLEX*16, INTENT(IN)    :: dmat_old(nat,4,-3:3,-3:3,3)
        REAL*8, INTENT(INOUT)     :: alx(nat,4,3),aly(nat,4,3),alz(nat,4,3)
        COMPLEX*16, INTENT(INOUT) :: dmat(nat,4,-3:3,-3:3,3)
        COMPLEX*16 ROKOLD,ROKNEW
        INTEGER LM(2,NCOM,NAT)
        DIMENSION CLMOLD(NRAD,NCOM,NAT,2),CLMNEW(NRAD,NCOM,NAT,2),       &
        ROKNEW(nq3,2),ROKOLD(nq3,2), dmixout(4)
        INTEGER mult(*)
        DIMENSION LMMAX(NAT),JRI(NAT),SCALE(*)
        integer  KVECVL(3,nq3)
        real*8 MSECINFO(20)
!
        CHARACTER*80 fname,fname2
        CHARACTER*5 MIX
        real*8 ZZ(*)
        real*8,allocatable :: PS(:,:),SB(:), Y(:,:),S(:,:)
        real*8,allocatable :: YY(:,:), FCURRENT(:),PCURRENT(:)
        real*8,allocatable :: SS(:,:),ProjHist(:)
        real*8,allocatable :: FHIST(:), PWHIST(:), CLMHIST(:),YHIST(:)
        logical inversion,reset,there, usererror
        character*4 MODUS, VERSION
        parameter (VERSION='VER3')
        parameter (ONE=1.D0, ZERO=0.0D0, dbase=0.005D0)
!
!       Debug/Full Version
#ifdef Full
        real*8,allocatable :: SY(:,:), Spart(:,:),Ypart(:,:)
        real*8,allocatable :: SY2(:,:), YY2(:,:), SS2(:,:)
        real*8,allocatable :: S2(:,:), Y2(:,:)
#endif
!
!       CHARGE DENSITY MIXING USING Various Methods
!       Authors, L. D. Marks and R. Luke
!       Version 4, November 2007
!
!       Update History
!       May 20th
!               Increased lower limit in Pratt to 0.01 -- otherwise below Wolf condition
!               Cleaned up a couple of format statements
!       May 21st
!               Lined up format 1002
!               Added back writing of :RED (sreduction, qmx)
!               Corrected reduction -- should not matter
!       May 24th
!               Gone to full internal control
!       November 2007
!               Added dynamic rescaling
!
        MIX='MSEC1'
!       Set Defaults
        DIAG =5D-4
        rtrap=0.1
        IDSCALE=0
        iscl=0
        icond=1
!       Trap too large/small values unless in Full mode
#ifndef Full
        if(qmx_input .ge. 0.6)then
                write(21,*)':WARNING: Mixing parameter may be too large and greedy'
                usererror=.true.
        else if(qmx_input .le. 0.025)then
                write(21,*)':WARNING: Mixing parameter may be too small'
                usererror=.true.
        endif
        if(nmax .lt. 1)then
                write(21,*)':WARNING: Number of memory steps too small, using 8'
                nmax=8
                usererror=.true.
        else if(nmax .lt. 4)then
                write(21,*)':WARNING: Number of memory steps many be too small'
                usererror=.true.
        endif
#endif
!       Read extra controls
        inquire(file='Tester', exist=there)
        if(there)then
                open(unit=1,file='Tester')
                read(1,*,err=983,end=983)IDTMP,ISCL, rtrap
                IDSCALE=IDTMP
                read(1,*,err=984,end=984)Tmp1
                DIAG = Tmp1
#ifdef Full
                read(1,*,err=984,end=984)itmp
                icond=itmp
#endif
                goto 984
983             IDSCALE=0
                rtrap =0.1
                ISCL=0
984             continue
        endif
!
!       Put in some traps
        DIAG=max(DIAG,1D-12)
        DIAG=min(DIAG,0.1D0)
        if(IDSCALE .lt. 0)IDSCALE=0
        if(IDSCALE .gt. 3)IDSCALE=0
!
!       Output the Options
        write(21,2001)' '
        write(21,2001)'******************************************************'
        write(21,2001)'* MULTISECANT MIXING OPTIONS                         *'
        write(21,2001)'* NFL Limit ',rtrap,'                               *'
        write(21,2001)'* Diag Parm ',diag,'                               *'
        if(idscale .eq. 0) &
        write(21,2001)'* Standard Increase/Decrease Bounds                  *'
        if(idscale .eq. 1) &
        write(21,2001)'* Standard Bounds with Old Wien Limit                *'
        if(idscale .eq. 2) &
        write(21,2001)'* Very Aggressive (Greedy) Bound                     *'
        if(idscale .eq. 3) &
        write(21,2001)'* Aggressive Bound with Old Wien Limit               *'
        if(icond .ne. 0) then
        write(21,2001)'* Plane Wave SQRT Dynamic Rescaling                  *'
        else
        write(21,2001)'* Full Plane Wave Dynamic Rescaling                  *'
        endif
         nmaxx=17
        if(nmax.gt.nmaxx) then
         write(21,'("* NMAX gt ",i2," not supported and automatically reduced *")') nmaxx
         nmax=nmaxx-1
        endif
        write(21,2002)'* Max Number of Memory Steps ',NMAX,'                    *'
       write(21,2001)'******************************************************'
        write(21,2001)' '
2002    format(a,i4,a)
!       Determine the total Size
!
        MAXMQ=NQ3
        idmat=0
        DO idm=1,ndm
                DO  JATOM=1,NAT
                        IF(jatom1(jatom,1,idm).EQ.0) EXIT
                        DO iorb=1,4
                                IF(jatom1(jatom,iorb,idm).EQ.0) EXIT
                                l=ll(jatom,iorb,idm)
                                idmat=idmat+3+(l*l+1)*(l*l+1)*2
                        ENDDO
                ENDDO
        ENDDO
        if(inversion)then
                NBSIZE=jspin*nq3+NRAD*SUM(LMMAX(1:NAT))*jspin+idmat
        else
                NBSIZE=jspin*nq3*2+NRAD*SUM(LMMAX(1:NAT))*jspin+idmat
        endif
        allocate ( PS(NBSIZE,2))
!
!       Scaling for Plane waves -- has an effect !!!!
        scl_plane=sqrt(VI)*scl1
!
!       Pack values into PS Array
        call PackValues(PS,ROKNEW,ROKOLD,CLMOLD,CLMNEW, &
             alx_old,aly_old,alz_old,alx,aly,alz,dmat,dmat_old, &
             NPLANE,JSPIN,NAT,JRI,LMMAX,JATOM,NDM,jatom1,ll,nq3,              &
             SPLANE, SCHARGE, scl_plane, MAXMIX,NBSIZE,MaxMQ,inversion, &
             TPLANE, TCharge)
!
!       Splane is the mean-square distance in real space
!       The factor of ndif puts it on a comparable scale to :DIS, :ENE
        if(.not.inversion)Splane=Splane*0.5D0
        if(Splane.gt.1d-20)Splane=sqrt(Splane)*ndif
        write(21,210)':PLANE:  INTERSTITIAL TOTAL ',Tplane, ' DISTAN ',Splane
        write(21,210)':CHARG:  CLM CHARGE   TOTAL ',TCharge,' DISTAN ',SCharge
!
!       Apply general bounds based upon residues and type of atoms
        call stepbound(sreduction,qmx_input,qmx,qtot,splane,scharge,dbase)
        write(21,4141)sreduction,qmx
!
!--------------------------------------------------------------------
!       Start Collecting History about previous steps
!       Arrange into S, Y matrices as we go
!--------------------------------------------------------------------
!
        RedOld=1
        RedPred=1
        read(31,err=6001,end=6001)dmix,niter,modus
        dmix_last=dmix
!       Get prior step and predicted reduction
        read(31,err=6002,end=6002)Tmp1
        RedOld=Tmp1
        goto 6002
6001    niter=0
6002    niter=niter+1
!       MEMALL is 1+the total number of prior iterations
        MEMALL=MIN(NMAX+1,NITER)
        ISKIP=NITER-MEMALL
        MEMORY=MEMALL-1
!       For Memory=1, DIAG should be zero (small)
        if(MEMORY .EQ. 1)DIAG=1D-12
!
        allocate (SB(MAXMIX))
        allocate (FCURRENT(MAXMIX), PCURRENT(MAXMIX))
!
!       Current point, not most efficient
        PCURRENT(1:MAXMIX)=PS(1:MAXMIX,1)
        FCURRENT(1:MAXMIX)=PS(1:MAXMIX,2)-PCURRENT(1:MAXMIX)
        deallocate (PS)
!
        allocate (FHIST(MEMALL), PWHIST(MEMALL), CLMHIST(MEMALL))
        allocate (YHIST(MEMALL) )
!
        FHIST(MEMALL)  =dot_product(FCURRENT,FCURRENT)
        PWHIST(MEMALL) =dot_product(FCURRENT(1:NPLANE),FCURRENT(1:NPLANE))
        CLMHIST(MEMALL)=FHIST(MEMALL)-PWHIST(MEMALL)
        if(MEMORY .gt. 0)then
                allocate (Y(MAXMIX,MEMORY), S(MAXMIX,MEMORY))
                allocate (YY(MEMORY,MEMORY),SS(MEMORY,MEMORY))
#ifdef Full
!      Full version
                allocate (SY(MEMORY,MEMORY))
                allocate (SY2(MEMORY,MEMORY),YY2(MEMORY,MEMORY), SS2(MEMORY,MEMORY))
                allocate (Y2(MAXMIX,MEMORY),S2(MAXMIX,MEMORY))
#endif
!PB
  close (32)
!                if(iskip.gt.0)write(21,*)'Skipping ',iskip,' Total ',niter
                DO N=1,MEMORY
                   call open_broydfile(fname,iskip+n,fname2)
!     print*, 'open ',fname2
                   open(32,file=fname2,status='old',form='unformatted')
                   read(32)SB
                   S(1:MAXMIX,N)=PCURRENT(1:MAXMIX)-SB(1:MAXMIX)
                   read(32)SB
                   close (32)
                   Y(1:MAXMIX,N)=FCURRENT(1:MAXMIX)-SB(1:MAXMIX)
                   FHIST(N)  =dot_product(SB(1:MAXMIX),SB(1:MAXMIX))
                   PWHIST(N) =dot_product(SB(1:NPLANE),SB(1:NPLANE))
                   CLMHIST(N)=FHIST(N)-PWHIST(N)
                enddo
        endif
!
#ifdef Full
!       Create the terms for conventional BB & GB
        DO N=1,MEMORY-1
                S2(1:MAXMIX,N)=S(1:MAXMIX,N+1)-S(1:MAXMIX,N)
                Y2(1:MAXMIX,N)=Y(1:MAXMIX,N+1)-Y(1:MAXMIX,N)
        ENDDO
        if(MEMORY .gt. 0)then
                S2(1:MAXMIX,MEMORY)=S(1:MAXMIX,MEMORY)
                Y2(1:MAXMIX,MEMORY)=Y(1:MAXMIX,MEMORY)
        endif
#endif
!
                   call open_broydfile(fname,iskip+memory+1,fname2)
                   open(32,file=fname2,status='unknown',form='unformatted')
                   write(32)PCURRENT
                   write(32)FCURRENT
                   close(32)
                   if(niter-nmaxx.ge.0) then
                      call open_broydfile(fname,niter-nmaxx+1,fname2)
                      OPEN (32,FILE=FNAME2,STATUS='REPLACE')
                      close (32,STATUS='DELETE')
                   endif
!
!--------------------------------------------------------------------
!       End Collecting History about previous steps
!--------------------------------------------------------------------
!
        if(MEMORY .eq. 0)goto 119
!
!       On-the-fly rescaling of PW versus CLM
!       Create a running average of the distances
        PWAVE=0.
        CLAVE=0.
        do N=1,MEMALL
                PWAVE=PWAVE+sqrt(PWHIST(N)/FHIST(N))
                CLAVE=CLAVE+sqrt(CLMHIST(N)/FHIST(N))
        enddo
!
!       For the PW's rescale so residue is comparable to that of CLMs
!       Default is icond=1, take sqrt
!       This makes sense because relative weights appear as Rescale**2 in algorithm
        Rescale = CLAVE/PWAVE
        if(icond .gt. 0)Rescale=sqrt(Rescale)
        MSECINFO(1)=Rescale
!
        write(21,1002)' Dynamic rescale ',rescale
        PCURRENT(1:NPLANE)=PCURRENT(1:NPLANE)*rescale
        FCURRENT(1:NPLANE)=FCURRENT(1:NPLANE)*rescale
        Y(1:NPLANE,1:MEMORY)=Y(1:NPLANE,1:MEMORY)*rescale
        S(1:NPLANE,1:MEMORY)=S(1:NPLANE,1:MEMORY)*rescale
#ifdef Full
        Y2(1:NPLANE,1:MEMORY)=Y2(1:NPLANE,1:MEMORY)*rescale
        S2(1:NPLANE,1:MEMORY)=S2(1:NPLANE,1:MEMORY)*rescale
#endif
        PWHIST(1:MEMALL)=PWHIST(1:MEMALL)*rescale*rescale
        FHIST(1:MEMALL)=PWHIST(1:MEMALL)+CLMHIST(1:MEMALL)
        scl_plane=scl_plane*rescale
!
!--------------------------------------------------------------------
!       Renormalize the S & Y vectors
!--------------------------------------------------------------------
!
        do J=1,MEMORY
           T1=0.
!
           DO K=1,MAXMIX
                T1=T1+Y(K,J)*Y(K,J)
           enddo
           YHIST(J) = T1
!       Renormalize
           T1=1.D0/sqrt(T1)
           DO K=1,MAXMIX
                S(K,J)=S(K,J)*T1
                Y(K,J)=Y(K,J)*T1
           ENDDO
#ifdef Full
!       Renormalize for BB, GB
           T1=dot_product(Y2(1:MAXMIX,J),Y2(1:MAXMIX,J))
           T1=1.D0/sqrt(T1)
           Y2(1:MAXMIX,J)=Y2(1:MAXMIX,J)
           S2(1:MAXMIX,J)=S2(1:MAXMIX,J)
#endif
!
        enddo
!
!--------------------------------------------------------------------
!       End Renormalizing the S & Y vectors
!--------------------------------------------------------------------
!
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
#ifdef Full
!       For full version
          DO K=1,MEMORY
                SS2(J,K)=dot_product(S2(1:MAXMIX,J),S2(1:MAXMIX,K))
                YY2(J,K)=dot_product(Y2(1:MAXMIX,J),Y2(1:MAXMIX,K))
                SY2(J,K)=dot_product(S2(1:MAXMIX,J),Y2(1:MAXMIX,K))
                SY(J,K) =dot_product(S (1:MAXMIX,J),Y (1:MAXMIX,K))
          ENDDO
#endif
!
        enddo
!       Do transpose part as well
        DO J=1,MEMORY
                DO K=1,J-1
                        SS(K,J)=SS(J,K)
                        YY(K,J)=YY(J,K)
                enddo
        enddo
!
#ifdef Full
!       Get the Condition Number of the matrices
!       Useful diagnostic, although these are horribly big!
!       CondYY =GetCondNo(YY ,MEMORY,DIAG)
!       write(21,1002)' Cond Number ',DIAG,CondYY
!       CondSS =GetCondNo(SS ,MEMORY,DIAG)
!       write(21,1002)' Cond Numbers ',CondYY,Condss
#endif
!
        MUSE=MEMORY
!
        write(21,1102) MUSE,iskip+MEMORY-MUSE
!
!--------------------------------------------------------------------
#ifdef Full
!       Bound the step here
        call LimitDMIX(Y,S,YY,FCURRENT,FHIST,YHIST,PWHIST,CLMHIST,MAXMIX,MEMALL,MEMORY, &
                             ascl1, qmx_input, dmix_last, qmx, dmixout, &
                             IDSCALE,Nplane, PM1, rtrap, Dbase, &
                             ISCL, RedPred, RedOld, DIAG, MUSE, &
                             S2,Y2,YY2,MSECINFO)
#else
        call LimitDMIX(Y,S,YY,FCURRENT,FHIST,YHIST,PWHIST,CLMHIST,MAXMIX,MEMALL,MEMORY, &
                             ascl1, qmx_input, dmix_last, qmx, dmixout, &
                             IDSCALE,Nplane, PM1, rtrap, Dbase, &
                             ISCL, RedPred, RedOld, DIAG, MUSE,MSECINFO)
#endif
!
        DMIXM=dmixout(1)
        DMIX=DMIXM
        qmx=dmixm
        if(nmax .lt. 10)then
        write(21,8001)MUSE,nmax,MSECINFO(1:4)
8001    format(':DIRM :  MEMORY ',i1,'/',i1,' RESCALE ',F6.3,' RED ',F6.3,' PRED ',F6.3,' NEXT ',F6.3)
        else
        write(21,8002)MUSE,nmax,MSECINFO(1:4)
8002    format(':DIRM :  MEMORY ',i2,'/',i2,' RESCALE ',F6.3,' RED ',F6.3,' PRED ',F6.3,' NEXT ',F6.3)
        endif
!
!--------------------------------------------------------------------
!
!       Is it the current version or not?
        if(MODUS .ne. VERSION) goto 119
!
!--------------------------------------------------------------------
!       Calculate the New Steps, currently with 3 modes
!
!--------------------------------------------------------------------
!       Set Diagonal Initial Matrix
!        call SetHZero(HZERO,DMIX,SCL1,MAXMIX,NPLANE)
!--------------------------------------------------------------------
!
!      Note: Y and S referenced to current position here -- works better
       IFAIL=0
#ifdef Full
       if(ISCL .eq. -2)then
        MIX ='BBROY'
        call  BBLuke(Y2,S2,YY2,FCURRENT,SB,MAXMIX,MEMORY,DMIXM)
       else if(ISCL .eq. -3)then
        MIX ='GBROY'
        call  GBLDM(Y2,S2,SY2,SS2,FCURRENT,SB,MAXMIX,MEMORY,DMIXM,IFAIL)
       else if(ISCL .eq. -4)then
        MIX ='MSEC2'
         call MSEC2(Y,S,SY,FCURRENT,SB,MAXMIX,MEMORY,DMIXM,IFAIL,DIAG)
       else if(ISCL .eq. -5)then
        MIX ='GRBO2'
        call  GBLDM(Y,S,SY,SS,FCURRENT,SB,MAXMIX,MEMORY,DMIXM,IFAIL)
       else
        call  MSEC1(Y,S,YY,FCURRENT,SB,MAXMIX,MEMORY,DMIXM,IFAIL,DIAG,MUSE)
       endif
#else
        call  MSEC1(Y,S,YY,FCURRENT,SB,MAXMIX,MEMORY,DMIXM,IFAIL,DIAG,MUSE)
#endif
!      Inversion bombed, go to Pratt
       if(IFAIL .ne. 0)then
                write(21,*)':WARNING: Inversion of Multi-Secant Matrix Failed'
                usererror=.true.
                goto 119
       endif
       write(21,1002)' DMIXM and Projections   ',DMIXM,PM1
!
!
        SB(1:MAXMIX)=PCURRENT(1:MAXMIX)+SB(1:MAXMIX)
        DMIXUSED=DMIXM
!
        DMIX=DMIXUSED
        QMX=DMIX
!
!--------------------------------------------------------------------
!       Done with calculating the updates
!
!       Consistency check
        nreduce=0
        reduction=1.D0
120     DOT=0
        BMOD=0
        PMOD=0
!
!       Calculate the angle between the Pratt & Broyden Steps
!       Also calculate the relative sizes
        DO K=1,NPLANE
                SPRATT=FCURRENT(K)*DMIX
                SBROY=(SB(K)-PCURRENT(K))
                DOT=SPRATT*SBROY+DOT
                BMOD=BMOD+SBROY*SBROY
                PMOD=PMOD+SPRATT*SPRATT
        enddo
!       Adjust to put plane waves on appropriate scale
        DOTP=DOT
        BMODP=BMOD
        PMODP=PMOD
        DTTP=DOT/sqrt(BMOD*PMOD)
        if(DTTP .gt. 1.0D0)then
                DTTP=0.0D0
        else if(DTTP .lt. -1.0D0) then
                DTTP=180.D0
        else
                DTTP=acos(DTTP)*90./asin(1.D0)
        endif
        DO K=1+NPLANE,MAXMIX
                SPRATT=FCURRENT(K)*DMIX
                SBROY=(SB(K)-PCURRENT(K))
                DOT=SPRATT*SBROY+DOT
                BMOD=BMOD+SBROY*SBROY
                PMOD=PMOD+SPRATT*SPRATT
        ENDDO
!       Traps based upon the relative sizes and the angle
!
!       Test the angle
!       If it is too big Broyden method may have collapsed
        BMOD=sqrt(BMOD)
        PMOD=sqrt(PMOD)
        DTT=DOT/(BMOD*PMOD)
        if(DTT .gt. 1.0D0)then
                DTT=0.0D0
        else if(DTT .lt. -1.0D0) then
                DTT=180.D0
        else
                DTT=acos(DTT)*90./asin(1.D0)
        endif
!
!       Trap too large a step
        BLIMIT=2.0*PMOD/dmix                  ! 2 times a full Pratt step
!       Lower bound for reduction
        RedLower = 0.25D0                     ! Could use 1/3, be conservative for production version
!
!       Ensure that we don't go too small
        BLIMIT=max(BLIMIT,5.0D-2)
!
!       Skip trap if broyden is very small
        if((bmod .lt. 1.D-2).or.(reduction .lt.(RedLower+1D-5)))goto 122
!
        if(BMOD .gt. BLIMIT)then
!                       Reduce the step
                        nreduce=nreduce+1
!                       In case something went wrong
                        if(nreduce .gt. 1)then
                                X1=0.95D0
                        else
                                X1=0.99999999D0*BLIMIT/BMOD
                        endif
!       Trap too severe a reduction
                        rednew=max(RedLower,reduction*X1)
                        X1=rednew/reduction
                        X2=(1.D0-X1)
                        reduction=rednew
                        do K=1,maxmix
                                SB(K)=SB(K)*X1+PCURRENT(K)*X2
                        enddo
                        goto 120
        endif
!
!       Done with traps, move to convert values back
122     continue
        write(21,982)sqrt(Bmodp),sqrt(pmodp),dttp
        if(reduction .lt. 0.9)then
            Write(21,981)BMOD,PMOD,DTT,reduction
        else
            Write(21,98) BMOD,PMOD,DTT
        endif

!--------------------------------------------------------------------
!       Put positions back into the arrays that they came from
!--------------------------------------------------------------------
!
100     CONTINUE
        Call UnPackValues(ROKNEW,CLMNEW, alx,aly,alz,dmat, &
             NPLANE,JSPIN,NAT,JRI,LMMAX,JATOM,NDM,jatom1,ll,nq3,              &
             MAXMIX,SB, scl_plane,MaxMQ,inversion)
        rewind 31
        write(31)DMIX,NITER,VERSION
        write(31)RedPred
        GOTO 121
!
!--------------------------------------------------------------------
!       PRATT MIXING FOR FIRST ITERATION, after angle gets too big
!       or something else goes wrong
!--------------------------------------------------------------------
!
119     CONTINUE
!       Note: We are resetting here, but this might not always be right
!       If we get here via a trap on the angle perhaps we should not reset
!
        niter=1
        call stepbound(sreduction,qmx_input,qmx,qtot,splane,scharge,dbase)
        DMIX=min(QMX,0.025D0)
        DMIX=max(DMIX,dbase)
        MIX='PRATT'
        QMX=DMIX
        DMIXUSED=DMIX
!
        DO K=1,MAXMIX
                SB(K)=PCURRENT(K)+DMIX*FCURRENT(K)
        ENDDO
!       Jump back to main output code
        GOTO 100

121     CONTINUE
!       Set QMX for output
        DMIX=DMIXUSED
        QMX =DMIX
!
!       Cleanup
        deallocate (FCURRENT, PCURRENT,SB)
        deallocate (FHIST, PWHIST, CLMHIST,YHIST)
        if(MEMORY .gt. 0)then
                deallocate (Y,S,YY,SS)
#ifdef Full
                deallocate (SY,SY2,YY2,,SS2,Y2,S2)
#endif
        endif
!
        RETURN
!
200     format(':PRATT STEP USED')
210 	FORMAT(A,F12.5,A,F11.7)
982     format(':DIRP :  |BROYD|=',D10.3,' |PRATT|=',D10.3, &
        ' ANGLE=',f6.1,' DEGREES')
981     format(':DIRB :  |BROYD|=',D10.3,' |PRATT|=',D10.3, &
        ' ANGLE=',f6.1,' DEGREES, STEP=',D10.3)
98      format(':DIRB :  |BROYD|=',D10.3,' |PRATT|=',D10.3, &
        ' ANGLE=',f6.1,' DEGREES')
1002    format(':INFO : ',a,10D11.3)
1102    format(':INFO :  Number of Memory Steps ',i4,' Skipping ',i4)
2001    format(a,F10.6,a)
2100    format(a,F10.6,2D14.6)
1022    FORMAT(10X,'DENSITY FOR ITERATION',I4,' PREPARED')
4141    format(':REDuction and DMIX in Broyd:',3f10.4,E14.5)
        END
