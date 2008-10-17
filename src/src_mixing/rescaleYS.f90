subroutine rescaleYS(n,S,Y)
use modmixermsec,only:noldstepsmax,PWHIST,FHIST,CLMHIST,yhist

implicit none
integer, intent(in)::n
real(8),intent(inout)::S(n,noldstepsmax),Y(n,noldstepsmax)
#ifdef none

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
#endif
 end subroutine
