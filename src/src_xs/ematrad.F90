!
!
!
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine ematrad (iq)
  use mod_qpoint, only: vql
  use m_writecmplxparts
  use m_getunit
      Use modmain
      Use modinput
      Use modxs
      Use m_getunit
      Implicit None
  ! arguments
      Integer, Intent (In) :: iq
  ! local variables
  integer :: un, un2, j
  character(256) :: fname, tmp1, tmpq1, tmpq2, tmpq3

      Integer :: is, ia, ias, nr, ir, igq
      Integer :: l1, l2, l3,lio,liomax
      Integer :: ilo, ilo1, ilo2, io, io1, io2
      Real (8) :: t1,ta,tb
      Integer :: lmax1, lmax2, lmax3
      Integer :: u11, u22, u33
  ! automatic arrays
      Real (8) :: r2 (nrmtmax), fr (nrmtmax), gr (nrmtmax), cf (3, &
     & nrmtmax)
  ! allocatable arrays
      Real (8), Allocatable :: jl (:, :), jhelp (:)
      integer,allocatable :: lio2l(:),lio2io(:)
!
      lmax1 = Max (input%xs%lmaxapwwf, lolmax)
      lmax2 = input%xs%lmaxemat
  ! lmax1 and lmax3 should be the same!
      lmax3 = lmax1
!
  ! allocate arrays for radial integrals and Bessel functions
      If (allocated(riaa)) deallocate (riaa)
      If (allocated(riloa)) deallocate (riloa)
      If (allocated(rilolo)) deallocate (rilolo)
      Allocate (riaa(0:lmax1, apwordmax, 0:lmax3, apwordmax, 0:lmax2, &
     & natmtot, ngq(iq)))
      Allocate (riloa(nlomax, 0:lmax3, apwordmax, 0:lmax2, natmtot, &
     & ngq(iq)))
      Allocate (rilolo(nlomax, nlomax, 0:lmax2, natmtot, ngq(iq)))
  ! allocate temporary arrays
!     Allocate (jl(0:lmax2, nrmtmax))
      Allocate (jl(nrmtmax,0:lmax2))
      Allocate (jhelp(0:lmax2))
      allocate(lio2l((lmax1+1)*apwordmax))
      allocate(lio2io((lmax1+1)*apwordmax))

      jl (:, :) = 0.d0
      jhelp (:) = 0.d0
  ! zero arrays for radial integrals
      riaa (:, :, :, :, :, :, :) = 0.d0
      riloa (:, :, :, :, :, :) = 0.d0
      rilolo (:, :, :, :, :) = 0.d0
!
      If (input%xs%dbglev .Gt. 1) Then
     ! APW-APW
         Call getunit (u11)
         Open (Unit=u11, File='IRADaa'//filext, Form='formatted', &
        & Action='write', Status='replace')
         Write (u11, '(a)') 'igq, ias, l1, io1, l3, io2, l2	 iraa'
         Write (u11, '(a)') '------------------------------------------&
        &-----------'
     ! lo-APW
         Call getunit (u22)
         Open (Unit=u22, File='IRADalo'//filext, Form='formatted', &
        & Action='write', Status='replace')
         Write (u22, '(a)') 'igq, ias, ilo, l1, l3, io, l2,	 iralo'
         Write (u22, '(a)') '------------------------------------------&
        &-----------'
     ! lo-lo
         Call getunit (u33)
         Open (Unit=u33, File='IRADlolo'//filext, Form='formatted', &
        & Action='write', Status='replace')
         Write (u33, '(a)') 'igq, ias, ilo1, l1, ilo2, l3, l2,   irlolo&
        &'
         Write (u33, '(a)') '------------------------------------------&
        &-----------'
      End If
if (.true.) then
!
  ! begin loop over G+q vectors
      Do igq = 1, ngq (iq)
     ! begin loop over species
         Do is = 1, nspecies
            nr = nrmt (is)
            Do ir = 1, nr
           ! calculate r^2
               r2 (ir) = spr (ir, is) ** 2
           ! calculate spherical Bessel functions of first kind j_l(|G+q|r_a)
               Call sbessel (lmax2, gqc(igq, iq)*spr(ir, is), jhelp)
               jl (ir,:) = jhelp (:)
            End Do

            lio=0
            Do l1 = 0, lmax1
              Do io1 = 1, apword (l1, is)
                lio=lio+1
                lio2l(lio)=l1
                lio2io(lio)=io1
              enddo
            enddo
            liomax=lio
        ! begin loop over atoms
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
           !----------------!
           !     APW-APW    !
           !----------------!
!               Do l2 = 0, lmax2
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(lio,l1,io1,l3,io2,ir,t1,fr,cf,gr)
!$OMP DO
#endif
               do lio=1,liomax
                 l1=lio2l(lio)
                 io1=lio2io(lio)
                     Do l3 = 0, lmax3
                        Do io2 = 1, apword (l3, is)
               Do l2 = 0, lmax2
                          Do ir = 1, nr
                            t1 = apwfr (ir, 1, io1, l1, ias) * apwfr (ir, 1, io2, l3, ias) * r2 (ir)
                            fr (ir) = t1 * jl (ir,l2)
                          End Do
                          Call fderiv (-1, nr, spr(1, is), fr, gr, cf)
                          riaa (l3, io2, l1, io1, l2, ias, igq) = gr (nr)
               End Do ! l2
                        End Do ! io2
                     End Do ! l3
!                  End Do ! io1
               End Do ! l1

#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL

#endif
!               End Do ! l2

!if( norm2( vql( :, iq) - (/0.25, 0.00, 0.00/)) .lt. input%structure%epslat) then
!  do l2 = 0, lmax3
!    do io2 = 1, apword( l2, is)
!      do l1 = 0, lmax1
!        do io1 = 1, apword( l1, is)
!          do l3 = 0, lmax2
!             write( *, '(2I3,3x,I3,3x,2I3,3x,2I3,3x,SP,F23.16)') is, ia, l3, l1, io1, l2, io2, riaa( l2, io2, l1, io1, l3, ias, igq)
!          end do
!        end do
!      end do
!    end do
!  end do
!end if
          
           !----------------------------!
           !     local-orbital-APW      !
           !----------------------------!
               Do ilo = 1, nlorb (is)
                  l1 = lorbl (ilo, is)
                  Do l3 = 0, lmax3
                     Do io = 1, apword (l3, is)
                        Do l2 = 0, lmax2
                           Do ir = 1, nr
                              t1 = lofr (ir, 1, ilo, ias) * apwfr (ir, 1, io, l3, ias) * r2 (ir)
                              fr (ir) = t1 * jl (ir,l2)
                           End Do
                           Call fderiv (-1, nr, spr(1, is), fr, gr, cf)
                           riloa (ilo, l3, io, l2, ias, igq) = gr (nr)
                        End Do ! l2
                     End Do ! io
                  End Do ! l3
               End Do ! ilo


           !------------------------------------!
           !     local-orbital-local-orbital    !
           !------------------------------------!
               Do ilo1 = 1, nlorb (is)
                  l1 = lorbl (ilo1, is)
                  Do ilo2 = 1, nlorb (is)
                     l3 = lorbl (ilo2, is)
                     Do l2 = 0, lmax2
                        Do ir = 1, nr
                           t1 = lofr (ir, 1, ilo1, ias) * lofr (ir, 1, ilo2, ias) * r2 (ir)
                           fr (ir) = t1 * jl (ir,l2)
                        End Do
                        Call fderiv (-1, nr, spr(1, is), fr, gr, cf)
                        rilolo (ilo1, ilo2, l2, ias, igq) = gr (nr)
                     End Do ! l2
                  End Do ! ilo2
               End Do ! ilo1
           
!      fname =''
!      tmp1 = ''
!      tmpq1 = ''
!      tmpq2 = ''
!      tmpq3 = ''
!      write(tmpq1, '(I4)') int(vql(1,iq)*1000)
!      write(tmpq2, '(I4)') int(vql(2,iq)*1000)
!      write(tmpq3, '(I4)') int(vql(3,iq)*1000)
!      write(tmp1,'(a,"_",a,"_",a)') trim(adjustl(tmpq1)),trim(adjustl(tmpq2)),trim(adjustl(tmpq3))
!      write(fname,'("ematri/",a,a,".OUT")') trim(adjustl('RILL')),trim(adjustl(tmp1))
!
!      call getunit(un)
!      open(unit=un, file=fname, action='write', status='replace')
!      call getunit(un2)
!      open(unit=un2, file=trim("sup"//fname), action='write', status='replace')
!  write(un2,*) "q vector"
!  write(un2, '(SP,E23.16)') vql(1,iq)
!  write(un2, '(SP,E23.16)') vql(2,iq)
!  write(un2, '(SP,E23.16)') vql(3,iq)
!  write(un2,*)
!  write(un2,*) "lmax2"
!  write(un2, '(I8)') lmax2
!  write(un2,*)
!  write(un2,'(a8,3(1x,a8))'), "ilo1","ilo2","l2"
!
!               Do ilo1 = 1, nlorb (is)
!                  Do ilo2 = 1, nlorb (is)
!                     Do l2 = 0, lmax2
!              write(un2, '(I8,3(1x,I8))') ilo1, ilo2, l2
!                     End Do ! l2
!                  End Do ! ilo2
!               End Do ! ilo1
!
!              j=0
!               Do ilo1 = 1, nlorb (is)
!                  Do ilo2 = 1, nlorb (is)
!                     Do l2 = 0, lmax2
!              j=j+1
!              write(un, '(SP,E23.16)') rilolo(ilo1, ilo2, l2, ias, igq)
!                     End Do ! l2
!                  End Do ! ilo2
!               End Do ! ilo1
!
!  write(un2,*)
!  write(un2,'("Number of entries: ",i10)'), j 
!
!  close(un2)
!      close(un)

!****************************************
! Debugging segment with output to files
!****************************************
               If (input%xs%dbglev .Gt. 1) Then
                   !----------------!
                   !     APW-APW    !
                   !----------------!
                       Do l1 = 0, lmax1
                          Do io1 = 1, apword (l1, is)
                             Do l3 = 0, lmax3
                                Do io2 = 1, apword (l3, is)
                                   Do l2 = 0, lmax2
                                      Write (u11, '(7i5, g18.10)') igq, ias, l1, io1, l3, io2, l2, riaa (l1, io1, l3, io2, l2, ias, igq)
                                   End Do
                                End Do ! io2
                             End Do ! l3
                          End Do ! io1
                       End Do ! l1
                   !----------------------------!
                   !     local-orbital-APW      !
                   !----------------------------!
                       Do ilo = 1, nlorb (is)
                          l1 = lorbl (ilo, is)
                          Do l3 = 0, lmax3
                             Do io = 1, apword (l3, is)
                                Do l2 = 0, lmax2
                                   Write (u22, '(7i5, g18.10)') igq, ias, ilo, l1, l3, io, l2, riloa (ilo, l3, io, l2, ias, igq) 
                                End Do ! l2
                             End Do ! io
                          End Do ! l3
                       End Do ! ilo
                   !------------------------------------!
                   !     local-orbital-local-orbital    !
                   !------------------------------------!
                       Do ilo1 = 1, nlorb (is)
                          l1 = lorbl (ilo1, is)
                          Do ilo2 = 1, nlorb (is)
                             l3 = lorbl (ilo2, is)
                             Do l2 = 0, lmax2
                                Write (u33, '(7i5, g18.10)') igq, ias, ilo1, l1, ilo2, l3, l2, rilolo (ilo1, ilo2, l2, ias, igq)
                             End Do ! l2
                          End Do ! ilo2
                       End Do ! ilo1
!**************************
! End of debugging segment
!**************************
                endif


           ! end loops over atoms and species
            End Do
         End Do
     ! end loop over G+q vectors
      End Do
endif
!
  ! deallocate
      Deallocate (jl, jhelp)
      deallocate(lio2l,lio2io)
      If (input%xs%dbglev .Gt. 1) Then
     ! close files
         Close (u11)
         Close (u22)
         Close (u33)
      End If
End Subroutine ematrad
