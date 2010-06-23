#ifdef Full
	Subroutine LimitDMIX(Y,S,YY,F,FHIST,YHIST,PWHIST,CLMHIST,MAXMIX,MEMALL,MEMORY, &
                             ascl1, qmx_input, dmix_last, qmx, dmixout, &
			     Mode,Nplane, PM1, rtrap, DBase, &
                             ISCL, RedPred, RedOld, DIAG,MUSE, &
                             S2,Y2,YY2,MSECINFO)
#else
        Subroutine LimitDMIX(Y,S,YY,F,FHIST,YHIST,PWHIST,CLMHIST,MAXMIX,MEMALL,MEMORY, &
                             ascl1, qmx_input, dmix_last, qmx, dmixout, &
                             Mode,Nplane, PM1, rtrap, DBase, &
                             ISCL, RedPred, RedOld, DIAG,MUSE,MSECINFO)
#endif
!
!       Bound the step based just upon prior history and NFL bound
!       Note: sections for MSGB, BB and GB are commented out unless -DFull is used during compilation
!
!       Modes
!               Mode            0       use increase/decrease bounds (default)
!                               1       add older Wien2k limit to 0
!                               2       greedy mode, use qmx_input unless diverging
!                               3       same as 2 with added older Wien2k limit
!       Mode 0 is the standard, 1 is slightly slower
!       Mode 2 can be better, but may also blowup so is not recommended. Mode 3 is similar.
!
	implicit real*8 (a-h,o-z)
	dimension FHIST(MEMORY),PWHIST(MEMORY),CLMHIST(MEMORY)
        dimension Y(MAXMIX,MEMORY),YY(MEMORY,MEMORY),F(MAXMIX),S(MAXMIX,MEMORY)
	dimension dmixout(4), YHIST(*)
        real*8 LowerBound, MSECINFO(*)
        parameter (LowerBound=0.025)
!       For Full version
#ifdef Full
        dimension S2(MAXMIX,MEMORY),Y2(MAXMIX,MEMORY)
        dimension YY2(MEMORY,MEMORY)
#endif
!
!--------------------------------------------------------------------
!       Limit based upon improvement
        PWLast   = PWHIST(MEMORY)
        CLMLast  = CLMHIST(MEMORY)
        PWThis   = PWHIST(MEMALL)
        CLMThis  = CLMHIST(MEMALL)
        RatioPW  = (PWThis*PWLast/(PWLast*PWLast+1D-50))
        RatioCLM = (CLMThis*CLMLast/(CLMLast*CLMLast+1D-50))
        RedGot   =  sqrt((RatioCLM + RatioPW)*0.5D0)
        Better   = sqrt(max(RatioPW, RatioCLM))
        AllBetter= sqrt(FHIST(MEMALL)/FHIST(MEMORY))
!
!       RedOld is what we predicted we would get as a reduction
!       RedGot is what we in fact achieved
!
!       Look at projected fractions for No Free Lunch Control
        PM1= Fprojmem (Y,S,F,YY,MAXMIX,MEMORY,PSTEP,DIAG,NPLANE,MUSE)
#ifdef Full
        if(ISCL .eq. -2)then
                PM1=BBProject(Y2,S2,YY2,F,MAXMIX,MEMORY)
        endif
#endif
        RedPred = PM1
!
#ifdef DEBUG
     write(*,221)RedGot,RedOld, RedPred,AllBetter
     format(':INFO :  Reduction ',F8.4,' Expected ',F8.4,' Next ',F8.4,' All ',F8.4)
#endif
        MSECINFO(2)=RedGot
        MSECINFO(3)=RedOld
        MSECINFO(4)=RedPred
        MSECINFO(5)=AllBetter
!
!	Limit based upon improvement
        qlimit1=qmx_input
        qlimit2=qmx_input
        if(mode .lt. 0)mode=0
        if(mode .gt. 3)mode=0
        if(mode .eq. 0)then
!              This appears to be the approximate best
                if(memory .eq. 1)dmix_last=max(dmix_last,0.015D0)
                if(AllBetter .lt. 1.0)then
                        qlimit1=min(2.0*dmix_last,dmix_last/AllBetter)
                else
                        qlimit1=max(0.5*dmix_last,dmix_last/AllBetter)
                endif
                if(memory .eq. 1) qlimit1=max(qlimit1,0.015D0)
                if(abs(qmx_input-qlimit1).lt.0.01)qlimit1=qmx_input
                qlimit2=qmx_input
        else if(mode .eq. 1)then
                if(memory .eq. 1)dmix_last=max(dmix_last,0.015D0)
                if(AllBetter .lt. 1.0)then
                        qlimit1=min(2.0*dmix_last,dmix_last/AllBetter)
                else
                        qlimit1=max(0.5*dmix_last,dmix_last/AllBetter)
                endif
                if(memory .eq. 1)qlimit1=max(qlimit1,0.015D0)
                if(abs(qmx_input-qlimit1).lt.0.01)qlimit1=qmx_input
!               Include Wien limit
                qlimit2=qmx
        else if (mode .eq. 2)then
                qlimit1=qmx_input
                qlimit2=qlimit1
        else if (mode .eq. 3)then
                qlimit1=qmx
                qlimit2=qmx
        endif
!
1000    continue
!
!       Combine the hard limits qlimit2
        dmix=min(qlimit1,qlimit2)
        if((qmx_input-dmix).lt.0.01)dmix=qmx_input
21      format(a,6D11.3)
!
        dmixout(1)=min(dmix,rtrap/max(PM1,0.001D0))
!
#ifdef DEBUG
        write(*,21)':INFO :  Bounds       ',qlimit1,qlimit2,rtrap/max(PM1,0.001D0),dmixout(1)
#endif
	do N=1,1
                t=max(dmixout(N),Dbase)
		dmixout(N)=t
	enddo
	return
	end
