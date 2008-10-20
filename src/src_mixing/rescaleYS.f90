subroutine rescaleYS(noldsteps,n,S,Y, potential,F)
use modmixermsec,only:noldstepsmax,PWHIST,FHIST,CLMHIST,yhist,FHIST,icond,scl_plane,MSECINFO
use modmain,only:ngrtot,lmmaxvr,nrmtmax,natmtot
implicit none
integer, intent(in)::n,noldsteps
real(8),intent(inout)::S(n,noldstepsmax),Y(n,noldstepsmax),potential(n),F(n)
real(8)::PWAVE,CLAVE,Rescale,T1
integer ::i,j,k,nmt
  nmt=lmmaxvr*nrmtmax*natmtot

    FHIST(noldsteps)  =dot_product(F,F)
    PWHIST(noldsteps) =dot_product(F(n-ngrtot:n),F(n-ngrtot:n))
    CLMHIST(noldsteps)=FHIST(noldsteps)-PWHIST(noldsteps)

		PWAVE=0.
        CLAVE=0.
        do i=1,noldsteps
                PWAVE=PWAVE+sqrt(PWHIST(i)/FHIST(i))
                CLAVE=CLAVE+sqrt(CLMHIST(i)/FHIST(i))
        enddo
!
!       For te PW's rescale so residue is comparable to that of CLMs
!       Default is icond=1, take sqrt
!       This makes sense because relative weights appear as Rescale**2 in algorithm
        Rescale = CLAVE/PWAVE
        if(icond .gt. 0)Rescale=sqrt(Rescale)
        MSECINFO(1)=Rescale
1002    format(':INFO : ',a,10D11.3)
        write(*,1002)' Dynamic rescale ',rescale
        potential(n-ngrtot:n)=potential(n-ngrtot:n)*rescale
        F(n-ngrtot:n)=F(n-ngrtot:n)*rescale
        Y(n-ngrtot:n,1:noldstepsmax)=Y(n-ngrtot:n,1:noldstepsmax)*rescale
        S(n-ngrtot:n,1:noldstepsmax)=S(n-ngrtot:n,1:noldstepsmax)*rescale
        PWHIST(1:noldsteps)=PWHIST(1:noldsteps)*rescale*rescale
        FHIST(1:noldsteps)=PWHIST(1:noldsteps)+CLMHIST(1:noldsteps)
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
