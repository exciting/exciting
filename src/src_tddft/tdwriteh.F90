
module m_tdwriteh
  implicit none
contains

  subroutine tdwriteh(un,iq)
    use modmain
    use modtddft
    implicit none
    ! arguments
    integer, intent(in) :: iq,un
    ! local variables
    character(10) :: dat,tim

    ! write prologue to file
    call date_and_time(date=dat,time=tim)
    write(unit1,*)
    write(unit1,'("## Date (YYYY-MM-DD): ",a4,"-",a2,"-",a2)') &
         dat(1:4),dat(5:6),dat(7:8)
    write(unit1,'("## Time (hh:mm:ss)  : ",a2,":",a2,":",a2)') &
         tim(1:2),tim(3:4),tim(5:6)
    write(unit1,'("# version           : ",i1.1,".",i1.1,".",i3.3)') version
    write(unit1,'("# build             : ",i5.5)') build_tddft
    write(unit1,'(a,3f12.6)') '# vql               :',vql(:,iq)
    write(unit1,'(a,3f12.6)') '# vqc               :',vqc(:,iq)
    write(unit1,'(a,2i8)') '# optcomp           :',optcomp(1,1),optcomp(2,1)
    write(unit1,'(a,i8)') '# fxctype           :',fxctype
    write(unit1,'(a,f12.6)') '# alphalrc          :',alphalrc
    write(unit1,'(a,f12.6)') '# alphalrcdyn       :',alphalrcdyn
    write(unit1,'(a,f12.6)') '# betalrcdyn        :',betalrcdyn
    write(unit1,'(a,l8)') '# aresdf            :',aresdf
    write(unit1,'(a,l8)') '# acont             :',acont
    write(unit1,'(a,i8)') '# nwacont           :',nwacont
    write(unit1,'(a,2f12.6)') '# brdtd (H,eV)      :',brdtd, h2ev*brdtd
    write(unit1,'(a,2f12.6)') '# scissor (H,eV)    :',scissor, h2ev*scissor
    write(unit1,'(a,i8)') '# nwdos             :',nwdos
    write(unit1,'(a,i8)') '# ngq               :',ngq(iq)
    write(unit1,'(a,f12.6)') '# gqmax             :',gqmax
    write(unit1,'(a,f12.6)') '# gmaxvr            :',gmaxvr
    write(unit1,'(a,f12.6)') '# rgkmax            :',rgkmax
    write(unit1,'(a,f12.6)') '# gkmax             :',gkmax
    write(unit1,'(a,3i8)') '# ngridk            :',ngridk
    write(unit1,'(a,3f12.6)') '# vkloff            :',vkloff
    write(unit1,'(a,l8)') '# reducek           :',reducek
    write(unit1,'(a,i8)') '# nmatmax           :',nmatmax
    write(unit1,'(a,i8)') '# ngkmax            :',ngkmax
    write(unit1,'(a,i8)') '# nlotot            :',nlotot
    write(unit1,'(a,i8)') '# nstval            :',nstval
    write(unit1,'(a,i8)') '# nstcon            :',nstcon
    write(unit1,'(a,i8)') '# nstsv             :',nstsv
    write(unit1,'(a,i8)') '# lmaxapw           :',lmaxapw
    write(unit1,'(a,i8)') '# lmaxapwtd         :',lmaxapwtd
    write(unit1,'(a,i8)') '# lmaxmat           :',lmaxmat
    write(unit1,'(a,i8)') '# lmaxvr            :',lmaxvr
    write(unit1,'(a,i8)') '# lmaxinr           :',lmaxinr
    write(unit1,'(a,i8)') '# lolmax            :',lolmax
    write(unit1,'(a,i8)') '# lmaxemat          :',lmaxemat
    write(unit1,'(a,i8)') '# lradstp           :',lradstp
    write(unit1,'(a,l8)') '# tevout            :',tevout
    write(unit1,'(a,l8)') '# nosym             :',nosym
    write(unit1,'(a,l8)') '# symwings          :',symwings
    write(unit1,'(a,l8)') '# tsymdfq0dn        :',tsymdfq0dn
    write(unit1,'(a,9f12.6)') '# symdfq0 (row1)    :',symdfq0(1,:)
    write(unit1,'(a,9f12.6)') '# symdfq0 (row2)    :',symdfq0(2,:)
    write(unit1,'(a,9f12.6)') '# symdfq0 (row3)    :',symdfq0(3,:)
    write(unit1,*)

  end subroutine tdwriteh

end module m_tdwriteh
