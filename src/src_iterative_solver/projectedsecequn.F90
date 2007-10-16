subroutine projectedsecequn(m,hp,op,evecp,evalp)
	
use modmain
implicit none
!update evecfv,eval  acording to sv
integer,intent(in)::m
complex(8),intent(in)::hp(2*m*(2*m+1)/2),op(2*m*(2*m+1)/2)
complex(8),intent(out)::evecp(2*m,m)
real(8),intent(out)::evalp(m)
complex(8)::work(2*2*m)
real(8):: abstol,rwork(7*2*m),vl,vu
integer nfound,iwork(5*2*m),ifail(2*m),info
real(8) dlamch
external dlamch

abstol=20.d0*dlamch('S')

!ZHPGVX(ITYPE, JOBZ, RANGE, UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, RWORK, IWORK, IFAIL, INFO)
call zhpgvx(1,'V','I','U',2*m,hp,op,vl,vu,1,m,abstol,nfound,evalp,evecp,2*m,work, &
 rwork,iwork,ifail,info)
	if(info.gt.0)then
		write(*,*) "ifail: ",ifail,"\ninfo: ",info,"\nnfound: ",nfound
 		stop
	end if

end  subroutine
