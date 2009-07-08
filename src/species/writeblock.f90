subroutine writeblock
use modspecies
use  mod_muffin_tin , only:epsedirac,epspotatom
implicit none
open(50,file=trim(spsymb)//trim(suffix)//'.in',action='WRITE',form='FORMATTED')
write(50,'(" ''",A,"''",T45,": spsymb")') trim(spsymb)
write(50,'(" ''",A,"''",T45,": spname")') trim(spname)
write(50,'(G14.6,T45,": spzn")') spzn
write(50,'(G18.10,T45,": spmass")') spmass
write(50,'(G14.6,2F10.4,I6,T45,": sprmin, rmt, sprmax, nrmt")') sprmin,rmt, &
 sprmax,nrmt
write(50,'(I4,T45,": spnst")') spnst
write(50,'(3I4,G14.6,L1,T45,": spn, spl, spk, spocc, spcore")') spn(1),spl(1), &
 spk(1),spocc(1),spcore(1)
do ist=2,spnst
  write(50,'(3I4,G14.6,L1)') spn(ist),spl(ist),spk(ist),spocc(ist),spcore(ist)
end do


! overall APW order
write(50,'(I4,T45,": apword")') apword
do io=1,apword
   write(50,'(F8.4,I4,"  ",L1,T45,": apwe0, apwdm, apwve")') boe,apwdm(io),apwve(io)
end do



if (apwordx.gt.0) then
  ! number of exceptions corresponds to number of l-values
  nlx=maxl+1
  write(50,'(I4,T45,": nlx")') nlx
  ! write the exceptions
  do l=0,maxl
    write(50,'(2I4,T45,": lorbl, lorbord")') l,apwordx
    do io=1,apwordx
      if (io.eq.1) then
        write(50,'(F8.4,I4,"  ",L1,T45,": apwe0, apwdm, apwve")') boe,apwdmx(io),apwvex(io)
      else
        write(50,'(F8.4,I4,"  ",L1,T45)') boe,apwdmx(io),apwvex(io)
      end if
    end do
  end do
else
  nlx=0
  write(50,'(I4,T45,": nlx")') nlx
end if

write(50,'(I4,T45,": nlorb")') nlorb

if (locorb) then
  ! write the local-orbitals
  do l=0,maxl
    write(50,'(2I4,T45,": lorbl, lorbord")') l,2
    if (searchlocorb) then
       write(50,'(F8.4,I4,"  ",L1,T45,": lorbe0, lorbdm, lorbve")') boe,0,.true.
       write(50,'(F8.4,I4,"  ",L1)') boe,1,.true.
    else
       write(50,'(F8.4,I4,"  ",L1,T45,": lorbe0, lorbdm, lorbve")') boe,0,.false.
       write(50,'(F8.4,I4,"  ",L1)') boe,1,.false.
    end if
  end do
end if

if (locorbsc) then
  do ist=1,spnst
    if (.not.spcore(ist)) then
      if ((spl(ist).eq.0).or.(spl(ist).eq.spk(ist))) then
        if (eval(ist).lt.esccut) then
          write(50,'(2I4,T45,": lorbl, lorbord")') spl(ist),3
          write(50,'(F8.4,I4,"  ",L1,T45,": lorbe0, lorbdm, lorbve")') boe,0, &
           .false.
          write(50,'(F8.4,I4,"  ",L1)') boe,1,.false.
          write(50,'(F8.4,I4,"  ",L1)') eval(ist)+0.5d0*boe,0,.true.
        end if
      end if
    end if
  end do
end if

write(50,*)
write(50,'("# Exciting code version : ",a)') version
write(50,'("# Description of method : ",a)') trim(apwdescr)

close(50)
end subroutine
