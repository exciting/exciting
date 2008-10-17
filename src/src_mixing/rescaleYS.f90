subroutine rescaleYS(MEMALL,n,S,Y, PCURRENT,FCURRENT)
use modmixermsec,only:noldstepsmax,PWHIST,FHIST,CLMHIST,yhist,FHIST,icond,scl_plane
use modmain,only:ngrtot
implicit none
integer, intent(in)::n,MEMALL
real(8),intent(inout)::S(n,noldstepsmax),Y(n,noldstepsmax),PCURRENT(n),FCURRENT(n)
real(8)::PWAVE,CLAVE,Rescale,T1
integer ::i,j,k


		PWAVE=0.
        CLAVE=0.
        do i=1,MEMALL
                PWAVE=PWAVE+sqrt(PWHIST(i)/FHIST(i))
                CLAVE=CLAVE+sqrt(CLMHIST(i)/FHIST(i))
        enddo
!
!       For te PW's rescale so residue is comparable to that of CLMs
!       Default is icond=1, take sqrt
!       This makes sense because relative weights appear as Rescale**2 in algorithm
        Rescale = CLAVE/PWAVE
        if(icond .gt. 0)Rescale=sqrt(Rescale)
       ! MSECINFO(1)=Rescale
1002    format(':INFO : ',a,10D11.3)
        write(*,1002)' Dynamic rescale ',rescale
        PCURRENT(1:ngrtot)=PCURRENT(ngrtot:n)*rescale
        FCURRENT(ngrtot:n)=FCURRENT(ngrtot:n)*rescale
        Y(ngrtot:n,1:noldstepsmax)=Y(ngrtot:n,1:noldstepsmax)*rescale
        S(ngrtot:n,1:noldstepsmax)=S(ngrtot:n,1:noldstepsmax)*rescale
        PWHIST(1:MEMALL)=PWHIST(1:MEMALL)*rescale*rescale
        FHIST(1:MEMALL)=PWHIST(1:MEMALL)+CLMHIST(1:MEMALL)
        scl_plane=scl_plane*rescale

  do J=1,noldstepsmax
           T1=0.
!
           DO K=1,n
                T1=T1+Y(K,J)*Y(K,J)
           enddo
           YHIST(J) = T1
!       Renormalize
           T1=1.D0/sqrt(T1)
           DO K=1,n
                S(K,J)=S(K,J)*T1
                Y(K,J)=Y(K,J)*T1
           ENDDO
end do
 end subroutine
