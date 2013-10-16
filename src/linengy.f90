!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: linengy
! !INTERFACE:
!
!
Subroutine linengy
! !USES:
      Use modinput
      Use modmain
      Use scl_xml_out_Module
! !DESCRIPTION:
!   Calculates the new linearisation energies for both the APW and local-orbital
!   radial functions. See the routine {\tt findband}.
!
! !REVISION HISTORY:
!   Created May 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ia, ja, ias, jas
      Integer :: l, ilo, io1, io2
! automatic arrays
      Logical :: done (natmmax)
      Real (8) :: vr (nrmtmax)
      logical :: tfnd
      character(256) :: linenetype, fname
!
      integer :: nr, nn
      real(8) :: t0, t1, e, dl
      real(8), allocatable :: p0 (:), p1 (:), q0 (:), q1 (:)
      character(1024) :: message

      linenetype = trim(input%groundstate%findlinentype)

!________________________________
!     switch off LE search

      if (linenetype.eq.'no_search') then
        do is = 1, nspecies
          do ia = 1, natoms (is)
            ias = idxas (ia, is)
! just copying default values for apw
            do l = 0, input%groundstate%lmaxapw
              do io1 = 1, apword(l,is)
                apwe(io1,l,ias) = apwe0(io1,l,is)
              end do
            end do
! and lo
            do ilo = 1, nlorb(is)
              do io1 = 1, lorbord(ilo,is)
                lorbe(io1,ilo,ias) = lorbe0(io1,ilo,is)
              end do
            end do
          end do
        end do
        return
      end if

!________________________________
!     l-charge based scheme

      if (linenetype.eq.'lcharge') then
          if ((iscl>1).and.(currentconvergence>0.0001)) call lchargelinene
          return
      end if
      
! begin loops over atoms and species
      Do is = 1, nspecies
         done (:) = .False.
         Do ia = 1, natoms (is)
            If ( .Not. done(ia)) Then
               ias = idxas (ia, is)
               vr (1:nrmt(is)) = veffmt (1, 1:nrmt(is), ias) * y00
!-----------------------!
!     APW functions     !
!-----------------------!
               Do l = 0, input%groundstate%lmaxapw
                  Do io1 = 1, apword (l, is)
                     apwe (io1, l, ias) = apwe0 (io1, l, is)
                     If (apwve(io1, l, is)) Then
! check if previous radial functions have same default energies
                        Do io2 = 1, io1-1
                           If (apwve(io2, l, is)) Then
                              If (Abs(apwe0(io1, l, is)-apwe0(io2, l, is)) .Lt. 1.d-4) Then
                                 apwe (io1, l, ias) = apwe (io2, l, ias)
                                 Go To 10
                              End If
                           End If
                        End Do
! find the band energy starting from default
                        Call findband (linenetype, &
                       &  l, 0, nrmt(is), &
                       &  spr(:, is), vr, input%groundstate%deband, input%groundstate%epsband, &
                       &  apwe(io1, l, ias),tfnd)
                        if (.not.tfnd) then
                          call warning('Warning(linengy):')
                          write(message,'(" Linearisation energy not found")')
                          call warning(message)
                          write(message,'(" for species ",I4," and atom ",I4)') is, ia
                          call warning(message)
                          write(message,'(" APW angular momentum ",I4)') l
                          call warning(message)
                          write(message,'(" order ",I4)') io1
                          call warning(message)
                          write(message,'(" and s.c. loop ",I5)') iscl
                          call warning(message)
                        end if
                     else
                       if (input%groundstate%fermilinengy) &
                         apwe(io1,l,ias)=efermi + input%groundstate%dlinengyfermi
                     End If
10                   Continue
                  End Do
               End Do
!---------------------------------!
!     local-orbital functions     !
!---------------------------------!
               Do ilo = 1, nlorb (is)
                  Do io1 = 1, lorbord (ilo, is)
                     lorbe (io1, ilo, ias) = lorbe0 (io1, ilo, is) 
                     If (lorbve(io1, ilo, is)) Then
                        l = lorbl (ilo, is)
! check if previous radial functions have same default energies
                        Do io2 = 1, io1 - 1
                           If (lorbve(io2, ilo, is)) Then
                              If (Abs(lorbe0(io1, ilo, is)-lorbe0(io2, ilo, is)) .Lt. 1.d-4) Then
                                 lorbe (io1, ilo, ias) = lorbe (io2, ilo, ias)
                                 Go To 20
                              End If
                           End If
                        End Do
                        Do io2 = 1, apword (l, is)
                           If (apwve(io2, l, is)) Then
                              If (Abs(lorbe0(io1, ilo, is)-apwe0(io2, l, is)) .Lt. 1.d-4) Then
                                 lorbe (io1, ilo, ias) = apwe (io2, l, ias)
                                 Go To 20
                              End If
                           End If
                        End Do
! find the band energy starting from default
                        Call findband (linenetype, &
                       &  l, 0, nrmt(is), &
                       &  spr(:, is), vr, input%groundstate%deband, input%groundstate%epsband, &
                       &  lorbe(io1, ilo, ias),tfnd)
                        if (.not.tfnd) then
                          call warning('Warning(linengy):')
                          write(message,'(" Linearisation energy not found")')
                          call warning(message)
                          write(message,'(" for species ",I4," and atom ",I4)') is, ia
                          call warning(message)
                          write(message,'(" local-orbital ",I4)') ilo
                          call warning(message)
                          write(message,'(" order ",I4)') io1
                          call warning(message)
                          write(message,'(" and s.c. loop",I5)') iscl
                          call warning(message)
                        end if
                     else
                       if (input%groundstate%fermilinengy) &
                         lorbe(io1,ilo,ias)=efermi + input%groundstate%dlinengyfermi
                     End If
20                   Continue
                  End Do
               End Do
               done (ia) = .True.
! copy to equivalent atoms
               Do ja = 1, natoms (is)
                  If (( .Not. done(ja)) .And. (eqatoms(ia, ja, is))) Then
                     jas = idxas (ja, is)
                     Do l = 0, input%groundstate%lmaxapw
                        Do io1 = 1, apword (l, is)
                           apwe (io1, l, jas) = apwe (io1, l, ias)
                        End Do
                     End Do
                     Do ilo = 1, nlorb (is)
                        Do io1 = 1, lorbord (ilo, is)
                           lorbe (io1, ilo, jas) = lorbe (io1, ilo, ias)
                        End Do
                     End Do
                     done (ja) = .True.
                  End If
               End Do
            End If
         End Do ! ia
!
! Visualize the logarithmic derivatives if required
!
         if ((trim(linenetype).eq.'logderiv').and.(tlast)) then
           nr = nrmt(is)
           allocate(p0(nr),p1(nr),q0(nr),q1(nr))
           do l = 0, 4
             write(fname,'("dl_l=",i1,"_is=",i1,".dat")') l, is
             open(777,file=fname,action='write')
             e = -10.d0
             do while (e .le. 10.d0)
               call rschroddme(0, l, 0, e, &
              &  nr, spr(:, is), vr, nn, p0, p1, q0, q1)
               t0 = p0(nr)
               t1 = p1(nr)
               if (dabs(t0)>1.0d-6) then
             	 ! shifted value of the logarithmic derivative
                 dl = spr(nr,is)*t1/t0+(l+1)
                 write(777,*) e, dl
               end if
               e = e + input%groundstate%deband
             end do ! e
             close(777)
           end do ! l
           deallocate(p0,p1,q0,q1)
         end if
      End Do ! is

      Return
End Subroutine
!EOC
