!BOP
!
!!ROUTINE: testmixfun
!
!!INTERFACE:

subroutine test_mixfun
      
!!DESCRIPTION:
!
!This subroutine perform the test of the mixed functions, calculates the
!product of two eigenvectors directly and as a linear combination of mixed
!basis functions and writes both to disk for plotting.      
!
!!USES:
    use modmain
    use modgw
    use mod_rpath

!!LOCAL VARIABLES:
    implicit none
    integer(4) :: ib1, ib2
    integer(4) :: ikp, jkp
    integer(4) :: ik, jk, iq
      
!!EXTERNAL ROUTINES:
    external calcmpwipw
    external eptest
    external plotevecprod
    external plotevecmix

!EOP
!BOC
    call boxmsg(6,'-','TESTMIXFUN')
    
    call init_rpath(rpath,input%gw%at1,input%gw%at2)
  
    ik = input%gw%iik
    jk = input%gw%jjk
            
    ikp = kset%ik2ikp(ik)
    jkp = kset%ik2ikp(jk)
      
    write(*,*) 'Parameters:'
    write(*,*) 'first k-point number (iik): ', input%gw%iik
    write(*,*) 'second k-point number (jjk): ', input%gw%jjk
    write(*,*) 'corresponding first irreducible k-point number (ikp): ', ikp
    write(*,*) 'corresponding second irreducible k-point number (jkp): ', jkp
    write(*,*) 'lower bound for band number (ibmin): ', input%gw%ibmin
    write(*,*) 'upper bound for band number (ibmax): ', input%gw%ibmax
    write(*,*) 'lower bound for band number (ibmin2): ', input%gw%ibmin2
    write(*,*) 'upper bound for band number (ibmax2): ', input%gw%ibmax2
    write(*,*) 'atom 1 (at1): ', input%gw%at1
    write(*,*) 'atom 2 (at2): ', input%gw%at2
    write(*,*)
    
    do iq = 1, kqset%nkpt
      if (kqset%kqid(ik,iq).eq.jk) exit
    enddo
    write(*,*) 'Corresponding Q-point: iq=', iq
    write(*,'(a,3f12.4)') ' vql=', kqset%vql(:,iq) 
    
    call eptest(ik,jk,iq)
    call diagsgi(iq)
    
    do ib1 = input%gw%ibmin, input%gw%ibmax
      do ib2 = input%gw%ibmin2, input%gw%ibmax2
    
        call plotevecprod(ik,jk,ib1,ib2,input%gw%at1,input%gw%at2)
    
        call plotevecmix(iq,ik,jk,ib1,ib2,input%gw%at1,input%gw%at2)
    
      enddo !ib2 
    enddo !ib1
    
    return
end subroutine
!EOC      
