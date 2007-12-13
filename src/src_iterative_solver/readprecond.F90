subroutine       readprecond(ik,n,X,w)
  use modmain
use modmpi	
  integer, intent(in)::n,ik
  complex(8), intent(out)::X(nmatmax,nmatmax)
   real(8),intent(out)::w(nmatmax)
  !local variables
  character(256) ::filetag
  character(256), external:: outfilenamestring
  integer recl,koffset
  inquire(iolength=recl)X,w
  filetag="PRECONDMATRIX"
  if (splittfile.or.(rank.eq.0)) then
     open(70,file=outfilenamestring(filetag,ik),action='READ', &
          form='UNFORMATTED',access='DIRECT',recl=recl)
     if (splittfile) then
        koffset=ik-firstk(procofk(ik))+1
     else
        koffset =ik
     endif
     read(70,rec=koffset)X,w
  else
     write(*,*)"Error"
     stop
  endif
end subroutine readprecond