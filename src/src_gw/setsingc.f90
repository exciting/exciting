subroutine setsingc()
        
      use modmain
      use modgw,  only: fgw,nqptnr,singc1,singc2,Gamma,gammapoint
      use modmpi, only: rank

      real(8) :: beta
      real(8) :: f1
      real(8) :: f2
      real(8) :: intf1   ! BZ integral of the auxiliary function 1      
      real(8) :: intf2   ! BZ integral of the auxiliary function 2      
      real(8) :: sumf1   ! Auxiliary functions for 
      real(8) :: sumf2   ! BZ integrals at singular
                         ! $\Gamma$ point
      integer :: iq
      
      beta=(omega/(6.0d0*pi*pi))**(1.0d0/3.0d0)
      intf1=omega/(4.0d0*pi*pi*beta)
      intf2=omega/(4.0d0*pi*pi)*sqrt(pi/beta)

      sumf1=0.0d0
      sumf2=0.0d0
      do iq = 1, nqptnr
        Gamma=gammapoint(iq)
        call genauxf(iq,beta,f1,f2)
        sumf1 = sumf1 + f1
        sumf2 = sumf2 + f2
      enddo  

      singc1=intf1-sumf1/dble(nqptnr)
      singc2=intf2-sumf2/dble(nqptnr)

      if (rank==0) then
        call linmsg(fgw,'-','SETSINGC evaluation of the singularity correction factors')
        write(fgw,100) beta
        write(fgw,101) intf1, intf2
        write(fgw,102) sumf1, sumf2
        write(fgw,103) singc1, singc2
      end if

  100 format('parameter beta=',f18.12,/,30x,'q^(-1)',12x,'q^(-2)')
  101 format('Numerical integration ',2f18.12)
  102 format('Analitic integration  ',2f18.12)
  103 format('Correction factor     ',2f18.12)
      return
end subroutine 
