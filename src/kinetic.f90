!
!
!
! Copyright (C) 2014 A. Gulans
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: kinetic
! !INTERFACE:
!
!
Subroutine KineticEnergy(ik,evecfv,apwalm,ngp,vgpc,igpig)
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Calculates the kinetic energy directly.
!
! !REVISION HISTORY:
!   Created December 2014 (A. Gulans)
!EOP
!BOC
      Implicit None
      integer, intent (in) :: ngp,ik
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, natmtot)
      Complex (8), Intent (in) :: evecfv (nmatmax, nstfv)
      Integer, Intent (In) :: igpig (ngkmax)
      Real (8), Intent (In) :: vgpc (3, ngkmax)

! local variables

      Integer :: is, ia, ias, nr, ir, if1,if3,inonz,ireset1,ireset3,maxi,i,ist
      Integer :: l1, l2, l3, m2, lm2, m1, m3, lm1, lm3
      Integer :: ilo, ilo1, ilo2, io, io1, io2, nalo1, maxnlo
      Real (8) :: t1,t2,angular
      real (8), allocatable :: t_aa(:,:,:),t_alo(:,:),t_lolo(:,:)
      Real (8) :: rmtable(nrmtmax),r2inv(nrmtmax)
      complex(8) :: zsum
      Real (8) :: r2 (nrmtmax), fr (nrmtmax), gr (nrmtmax), cf (3, &
     & nrmtmax),a,rm,alpha
      parameter (alpha=1d0 / 137.03599911d0)
      complex(8), allocatable :: zm(:,:),zvec(:),zveclo(:),zfft(:),zfft2(:)
      complex(8) :: zdotc
      external zdotc
      logical :: applyiora

      integer :: ifg,ix,igk,l,m,lm,LOoffset
      real(8) :: Eiora
      
      engyknst(:,ik) = 0.d0

! Initialisation of some variables that exist just for the sake of convenience    

      if (input%groundstate%ValenceRelativity.ne.'none') then
        a=0.5d0*alpha**2
      else
        a=0d0
      endif
      applyiora=((input%groundstate%ValenceRelativity.eq.'iora*').or.(input%groundstate%ValenceRelativity.eq.'iora*'))
      if (applyiora) then
        write(*,*) 'Direct kinetic energy calculations with IORA are not implemented yet... '
        stop
! evalfv needs to be incorporated
      endif

! MT part

! APW-APW storage initialisation
      
      if (nlomax.ne.0) then
        Allocate (t_lolo( nlomax,nlomax))     
        Allocate (t_alo( apwordmax, nlomax))     
!        write(*,*) 'allocated'
        t_lolo=0d0
        t_alo=0d0
      endif
!      write(*,*) nlomax
!      read(*,*)
      Allocate (t_aa ( apwordmax, apwordmax, 0:input%groundstate%lmaxmat))
      t_aa=0d0

      allocate(zm(apwordmax*lmmaxapw,nstfv))
      allocate(zvec(apwordmax*lmmaxapw))
     
! begin loops over atoms and species
      LOoffset=ngp
      Do is = 1, nspecies
         
         nr = nrmt (is)
         Do ir = 1, nr
            r2 (ir) = spr (ir, is) ** 2
            r2inv(ir)=1d0/r2(ir)
         End Do
         if (haloijSize(is).ne.0) allocate(zveclo(haloijSize(is)))
         Do ia = 1, natoms (is)
           ias = idxas (ia, is)
           Do ir = 1, nr
             rmtable (ir) = 1d0/(1d0-0.5d0*alpha**2*veffmt (1, ir, ias)*y00)
           End Do
!---------------------------!
!     APW-APW integrals     !
!---------------------------!
            Do l1 = 0, input%groundstate%lmaxmat
              Do io1 = 1, apword (l1, is)
                Do io2 = 1, apword (l1, is)
                  angular=dble(l1*(l1+1))
                  if  (input%groundstate%ValenceRelativity.eq.'none') then
! non-relativistic part
                    Do ir = 1, nr
                      t1=apwfr(ir, 1, io1, l1, ias)*apwfr(ir, 1, io2, l1, ias)
                      t2=apwfr(ir, 2, io1, l1, ias)*apwfr(ir, 2, io2, l1, ias)
                      fr (ir) = r2(ir)*(0.5d0*t2 + 0.5d0*angular*t1*r2inv(ir))
                    End Do
                    Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                    t_aa ( io2, io1, l1)= gr (nr)   !*4d0*pi 
                  else
! zora
                    Do ir = 1, nr
                      t1=apwfr(ir, 1, io1, l1, ias)*apwfr(ir, 1, io2, l1, ias)
                      t2=apwfr(ir, 2, io1, l1, ias)*apwfr(ir, 2, io2, l1, ias)
                      fr (ir) = (0.5d0*t2*rmtable(ir) + 0.5d0*angular*t1*rmtable(ir)*r2inv(ir))*r2 (ir)
                    End Do
			    Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                    t_aa ( io2, io1, l1)= gr (nr) !*4d0*pi
                    if (applyiora) then
! iora(1) correction
                      Do ir = 1, nr
                        t1=apwfr(ir, 1, io1, l1, ias)*apwfr(ir, 1, io2, l1, ias)
                        t2=apwfr(ir, 2, io1, l1, ias)*apwfr(ir, 2, io2, l1, ias)
                        fr (ir) = (0.5d0*t2*rmtable(ir)**2 + 0.5d0*angular*t1*rmtable(ir)**2*r2inv(ir))*r2 (ir)
                      End Do
                      Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                      t_aa ( io2, io1, l1)= t_aa ( io2, io1, l1)-0.25d0*a*gr (nr) !-0.5d0*alpha**2*gr (nr)
                    endif
                  endif
                End Do
              End Do
            End Do

           call zgemm('T', &           ! TRANSA = 'C'  op( A ) = A**H.
                      'N', &           ! TRANSB = 'N'  op( B ) = B.
                      apwordmax*lmmaxapw, &          ! M ... rows of op( A ) = rows of C
                      nstfv, &           ! N ... cols of op( B ) = cols of C
                      ngp, &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      apwalm(1,1,1,ias), &           ! A
                      ngkmax,&           ! LDA ... leading dimension of A
                      evecfv, &           ! B
                      nmatmax, &          ! LDB ... leading dimension of B
                      zzero, &          ! beta
                      zm, &  ! C
                      apwordmax*lmmaxapw &      ! LDC ... leading dimension of C
                     )

              
!--------------------------------------!
!     local-orbital-APW integrals      !
!--------------------------------------!
            Do ilo = 1, nlorb (is)
               l1 = lorbl (ilo, is)
               Do l3 = 0, input%groundstate%lmaxmat
                  Do io = 1, apword (l3, is)
                     If (l1 .Eq. l3) Then
                        angular=dble(l1*(l1+1))
                        Do ir = 1, nr
                           t1=apwfr(ir, 1, io, l1, ias)*lofr(ir, 1, ilo, ias)
                           t2=apwfr(ir, 2, io, l1, ias)*lofr(ir, 2, ilo, ias)
                           fr (ir) = (0.5d0*t2*rmtable(ir) + 0.5d0*angular*t1*rmtable(ir)*r2inv(ir))*r2 (ir)
                        End Do
                        Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                        t_alo(io,ilo)=gr(nr)
                        if (applyiora) then
! iora(1) correction
                          Do ir = 1, nr
                            t1=apwfr(ir, 1, io, l1, ias)*lofr(ir, 1, ilo, ias)
                            t2=apwfr(ir, 2, io, l1, ias)*lofr(ir, 2, ilo, ias)
                            fr (ir) = (0.5d0*t2*rmtable(ir)**2 + 0.5d0*angular*t1*rmtable(ir)**2*r2inv(ir))*r2 (ir)
                          End Do
                          Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                          t_alo (io, ilo)= t_alo (io, ilo)-0.25d0*a*gr (nr) !-0.5d0*alpha**2*gr (nr)
                        endif
                     End If
                  End Do
               End Do
            End Do
!-----------------------------------------------!
!     local-orbital-local-orbital integrals     !
!-----------------------------------------------!
            Do ilo1 = 1, nlorb (is)
               l1 = lorbl (ilo1, is)
               Do ilo2 = 1, nlorb (is)
                  l3 = lorbl (ilo2, is)
                  If (l1 .Eq. l3) Then
                     angular=dble(l1*(l1+1))
                     Do ir = 1, nr
                        t1=lofr(ir, 1, ilo1, ias)*lofr(ir, 1, ilo2, ias)
                        t2=lofr(ir, 2, ilo1, ias)*lofr(ir, 2, ilo2, ias)
                        fr (ir) = (0.5d0*t2*rmtable(ir) + 0.5d0*angular*t1*rmtable(ir)*r2inv(ir))*r2 (ir)
                     End Do
                     Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                     t_lolo(ilo1,ilo2)=gr(nr)
                     if (applyiora) then
! iora(1) correction
                       Do ir = 1, nr
                         t1=lofr(ir, 1, ilo1, ias)*lofr(ir, 1, ilo2, ias)
                         t2=lofr(ir, 2, ilo1, ias)*lofr(ir, 2, ilo2, ias)
                         fr (ir) = (0.5d0*t2*rmtable(ir)**2 + 0.5d0*angular*t1*rmtable(ir)**2*r2inv(ir))*r2 (ir)
                       End Do
                       Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                       t_lolo (ilo1, ilo2)= t_lolo (ilo1, ilo2)-0.25d0*a*gr (nr) !-0.5d0*alpha**2*gr (nr)
                     endif 
                  End If
               End Do
            End Do

!        zax2=0d0
!       if3=0


        
do ist=1,nstfv
              zvec=zzero
              if (haloijSize(is).ne.0) zveclo=zzero
              Do l3 = 0, input%groundstate%lmaxmat
                Do m3 = - l3, l3
                  lm3 = idxlm (l3, m3)
                  Do io2 = 1, apword (l3, is)
                    Do io1 = 1, apword (l3, is)
                      zvec(apwordmax*(lm3-1)+io2)=zvec(apwordmax*(lm3-1)+io2)+t_aa(io1,io2,l3)*zm(apwordmax*(lm3-1)+io1,ist)
                    enddo
                  End Do
                End Do
              End Do

        if3=0
        do ilo = 1, nlorb (is)
          l=lorbl (ilo, is)
! LO-APW and APW-LO
if (.true.) then
          do m=-l,l
            lm=idxlm (l, m)
!            if3=if3+1
            do io1=1,apword(l,is)
              zvec(apwordmax*(lm-1)+io1)=zvec(apwordmax*(lm-1)+io1)+t_alo(io1,ilo)*evecfv(LOoffset+if3+m+l+1,ist)
              zveclo(if3+m+l+1)=zveclo(if3+m+l+1)+t_alo(io1,ilo)*zm(apwordmax*(lm-1)+io1,ist)
            enddo
          enddo
endif
! LO-LO
if (.true.) then
         if1=0
         Do ilo2 = 1, nlorb (is)
           If (lorbl(ilo2, is) .Eq. l) Then
             do m=-l,l
                zveclo(if3+m+l+1)=zveclo(if3+m+l+1)+t_lolo(ilo,ilo2)*evecfv(LOoffset+if1+m+l+1,ist)
              enddo
            endif
            if1=if1+2*lorbl(ilo2, is)+1
          enddo
endif
          if3=if3+2*l+1
        enddo
!        write(*,*) if3,haloijSize(is)
!        write(*,*) LOoffset
!        write(*,*) 
        engyknst(ist,ik)=engyknst(ist,ik)+zdotc(apwordmax*lmmaxapw,zm(1,ist),1,zvec,1)
        if (nlorb(is).ne.0) engyknst(ist,ik)=engyknst(ist,ik)+zdotc(haloijSize(is),evecfv(LOoffset+1,ist),1,zveclo,1)
!        write(*,*) zdotc(apwordmax*lmmaxapw,zm(1,ist),1,zvec,1)+zdotc(haloijSize(is),evecfv(LOoffset+1,ist),1,zveclo,1)
enddo


! end loops over atoms and species
         LOoffset=LOoffset+haloijSize(is)
         End Do
         if (haloijSize(is).ne.0) deallocate(zveclo)
      End Do
! cleaning up 
      deallocate(t_aa)
      if (nlomax.ne.0) deallocate(t_alo,t_lolo)
      deallocate(zm)
      deallocate(zvec)

      allocate(zfft(ngrtot),zfft2(ngrtot))
      allocate(zvec(ngp))
      do ist=1,nstfv
        do ix=1,3

          zfft=zzero
          Do igk = 1, ngp
            ifg = igfft (igpig(igk))
            zfft (ifg) = evecfv (igk, ist) *vgpc(ix, igk)
          End Do
          Call zfftifc (3, ngrid, 1, zfft)
if  (input%groundstate%ValenceRelativity.eq.'none') then
          Do ir = 1, ngrtot
            zfft (ir)=zfft (ir)*cfunir(ir)
          End Do
else
          if (applyiora) zfft2=zfft
          Do ir = 1, ngrtot
            zfft (ir)=zfft (ir)*cfunir(ir)/(1d0-0.5d0*alpha*alpha*veffir(ir))
          End Do
 
endif
          Call zfftifc (3, ngrid,-1, zfft)
          zvec=zzero
          do igk=1,ngp
            zvec(igk)=zfft(igfft(igpig(igk)))*vgpc(ix, igk)
          enddo
          engyknst(ist,ik)=engyknst(ist,ik)+0.5d0*zdotc(ngp,evecfv(1,ist),1,zvec,1)

if (applyiora) then
          zfft=zfft2
          Do ir = 1, ngrtot
            zfft (ir)=zfft (ir)*cfunir(ir)/(1d0-0.5d0*alpha*alpha*veffir(ir))**2
          End Do
          Call zfftifc (3, ngrid,-1, zfft)
          zvec=zzero
          do igk=1,ngp
            zvec(igk)=zfft(igfft(igpig(igk)))*vgpc(ix, igk)
          enddo
          engyknst(ist,ik)=engyknst(ist,ik)-0.25d0*a*zdotc(ngp,evecfv(1,ist),1,zvec,1) ! evalfv is missing
endif 

        enddo
      enddo

      deallocate(zvec)
      deallocate(zfft,zfft2)

      Return
End Subroutine
!EOC

