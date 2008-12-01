subroutine rescaleYS(noldsteps,n,S,Y, potential,residual)
use modmixermsec,only:noldstepsmax,PWHIST,FHIST,CLMHIST,yhist,FHIST,icond,scl_plane,MSECINFO
use mod_Gvector,only:ngrtot
use mod_muffin_tin,only:lmmaxvr,nrmtmax
use mod_atoms,only:natmtot
implicit none
integer, intent(in)::n,noldsteps
real(8),intent(inout)::S(n,noldstepsmax),Y(n,noldstepsmax)
real(8),intent(inout)::potential(n),residual(n)
real(8)::PWAVE,CLAVE,Rescale,T1
integer ::i,j,k,nmt,firstpw,lastpw

  nmt=lmmaxvr*nrmtmax*natmtot
firstpw=n-ngrtot
lastpw=n
    FHIST(noldsteps)  =dot_product(residual,residual)

    PWHIST(noldsteps) =dot_product(residual(firstpw:lastpw),residual(firstpw:lastpw))
    CLMHIST(noldsteps)=FHIST(noldsteps)-PWHIST(noldsteps)
!Preconditioner Omega_n Pg 21
#ifdef Full_Understanding_how_this_should_work
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
        !
        !rescale doesnt work yet in exciting
        Rescale=1
        !
        MSECINFO(1)=Rescale
1002    format(':INFO : ',a,10D11.3)
        write(*,1002)' Dynamic rescale ',rescale
        potential(firstpw:lastpw)=potential(firstpw:lastpw)*rescale
        residual(firstpw:lastpw)=residual(firstpw:lastpw)*rescale
        Y(firstpw:lastpw,1:noldstepsmax)=Y(firstpw:lastpw,1:noldstepsmax)*rescale
        S(firstpw:lastpw,1:noldstepsmax)=S(firstpw:lastpw,1:noldstepsmax)*rescale
        PWHIST(1:noldsteps)=PWHIST(1:noldsteps)*rescale*rescale
        FHIST(1:noldsteps)=PWHIST(1:noldsteps)+CLMHIST(1:noldsteps)
        scl_plane=scl_plane*rescale


#endif

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
