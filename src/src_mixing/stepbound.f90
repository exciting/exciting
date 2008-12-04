        subroutine stepbound(reduction)
        use modmixermsec,only:SCHARGE,splane,dbase,qmx_input,qmx,qtot
        implicit none
        real(8),intent(out)::reduction
        real(8)::limit,DSlope,PFACT
!
!       Simpler form
!       Set the limiting term based upon the maximum of
!               qtot:           The charge difference
!               splane:         The PW difference
!               Scharge:        The CLM difference
!               dbase:          Lower Bound
!
        parameter (DSlope =2.0D0)       ! How much to reduce exponentially
        parameter (PFACT  =3.5D0  )     ! Controls reduction in terms of limit
        limit=DSlope*max(qtot,splane/PFACT)
        reduction=0.1+exp(-limit)
        qmx=qmx_input*reduction
        if(qmx .lt. DBASE) qmx=DBASE
        qmx=min(qmx,qmx_input,1.0D0)

        return
        end

