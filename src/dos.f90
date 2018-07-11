! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: dos
! !INTERFACE:
Subroutine dos
! !USES:
      Use modinput
      Use modmain
      Use FoX_wxml
      Use modmpi, Only: rank, splittfile
! !DESCRIPTION:
!   Produces a total and partial density of states (DOS) for plotting. The total
!   DOS is written to the file {\tt TDOS.OUT} while the partial DOS is written
!   to the file {\tt PDOS\_Sss\_Aaaaa.OUT} for atom {\tt aaaa} of species
!   {\tt ss}. In the case of the partial DOS, each symmetrised
!   $(l,m)$-projection is written consecutively and separated by blank lines.
!   If the global variable {\tt lmirep} is {\tt .true.}, then the density matrix
!   from which the $(l,m)$-projections are obtained is first rotated into a
!   irreducible representation basis, i.e. one that block diagonalises all the
!   site symmetry matrices in the $Y_{lm}$ basis. Eigenvalues of a quasi-random
!   matrix in the $Y_{lm}$ basis which has been symmetrised with the site
!   symmetries are written to {\tt ELMIREP.OUT}. This allows for identification
!   of the irreducible representations of the site symmetries, for example $e_g$
!   or $t_{2g}$, by the degeneracies of the eigenvalues. In the plot, spin-up is
!   made positive and spin-down negative. See the routines {\tt gendmat} and
!   {\tt brzint}.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Logical :: tsqaz
      Integer :: lmax, lmmax, l, m, lm
      Integer :: ispn, jspn, is, ia, ias
      Integer :: ik, nsk (3), ist, iw, jst, n, i
      Real (8) :: dw, th, t1
      Real (8) :: v1 (3), v2 (3), v3 (3)
      Character (512) :: buffer
      Character(256) :: string
      Type (xmlf_t), Save :: xf
      Complex (8) su2 (2, 2), dm1 (2, 2), dm2 (2, 2)
      Character (256) :: fname
! allocatable arrays
      Real (8), Allocatable :: e (:, :, :)
      Real (8), Allocatable :: f (:, :)
      Real (8), Allocatable :: w (:)
      Real (8), Allocatable :: g(:,:), gp(:)
      Real (8), Allocatable :: edif(:,:,:), ej(:,:), fj(:,:)
! low precision for band character array saves memory
      Real (4), Allocatable :: bc (:, :, :, :, :)
      Real (8), Allocatable :: elm (:, :)
      Complex (8), Allocatable :: ulm (:, :, :)
      Complex (8), Allocatable :: a (:, :)
      Complex (8), Allocatable :: dmat (:, :, :, :, :)
      Complex (8), Allocatable :: sdmat (:, :, :, :)
      Complex (8), Allocatable :: apwalm (:, :, :, :, :)
      Complex (8), Allocatable :: evecfv (:, :, :)
      Complex (8), Allocatable :: evecsv (:, :)
      
      if (associated(input%groundstate%Hybrid)) then
        if (input%groundstate%Hybrid%exchangetypenumber == 1) then
          input%groundstate%stypenumber = -1
        end if
      end if

! initialise universal variables
      splittfile=.false.
      Call init0
      Call init1
      lmax = Min (4, input%groundstate%lmaxapw)
      lmmax = (lmax+1) ** 2
! allocate local arrays
      Allocate (e(nstsv, nkpt, nspinor))
      Allocate (f(nstsv, nkpt))
      Allocate (w(input%properties%dos%nwdos))
      Allocate (g(input%properties%dos%nwdos, nspinor))
      Allocate (gp(input%properties%dos%nwdos))
      Allocate (bc(lmmax, nspinor, natmtot, nstsv, nkpt))
      If (input%properties%dos%lmirep) Then
         Allocate (elm(lmmax, natmtot))
         Allocate (ulm(lmmax, lmmax, natmtot))
         Allocate (a(lmmax, lmmax))
      End If
      Allocate (dmat(lmmax, lmmax, nspinor, nspinor, nstsv))
      Allocate (sdmat(nspinor, nspinor, nstsv, nkpt))
      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv))
      Allocate (evecfv(nmatmax, nstfv, nspnfv))
      Allocate (evecsv(nstsv, nstsv))
! read density and potentials from file
        If (associated(input%groundstate%Hybrid)) Then
           If (input%groundstate%Hybrid%exchangetypenumber == 1) Then
! in case of HF hybrids use PBE potential
            string=filext
            filext='_PBE.OUT'
            Call readstate
            filext=string
           Else
               Call readstate
           End If
        Else         
           Call readstate
        End If 
! read Fermi energy from file
      Call readfermi
! find the new linearisation energies
      Call linengy
! generate the APW radial functions
      Call genapwfr
! generate the local-orbital radial functions
      Call genlofr
! update potential in case if HF Hybrids
        If (associated(input%groundstate%Hybrid)) Then
           If (input%groundstate%Hybrid%exchangetypenumber == 1) Then
               Call readstate
           End If
        End If 
! generate unitary matrices which convert the (l,m) basis into the irreducible
! representation basis of the symmetry group at each atomic site
      If (input%properties%dos%lmirep) Then
         Call genlmirep (lmax, lmmax, elm, ulm)
      End If
! compute the SU(2) operator used for rotating the density matrix to the
! desired spin-quantisation axis
      v1 (:) = input%properties%dos%sqados
      t1 = Sqrt (v1(1)**2+v1(2)**2+v1(3)**2)
      If (t1 .Le. input%structure%epslat) Then
        if (rank==0) then
          Write (*,*)
          Write (*, '("Error(dos): spin-quantisation axis (sqados) has zero length")')
          Write (*,*)
          Stop
        end if
      End If
      v1 (:) = v1 (:) / t1
      If (v1(3) .Ge. 1.d0-input%structure%epslat) Then
         tsqaz = .True.
      Else
         tsqaz = .False.
         v2 (1:2) = 0.d0
         v2 (3) = 1.d0
         Call r3cross (v1, v2, v3)
! note that the spin-quantisation axis is rotated, so the density matrix should
! be rotated in the opposite direction
         th = - Acos (v1(3))
         Call axangsu2 (v3, th, su2)
      End If

! loop over k-points
      Do ik = 1, nkpt
! get the eigenvalues/vectors from file
         Call getevalsv (vkl(1, ik), evalsv(1, ik))
         
         If (input%properties%dos%lmirep) Then
         
         Call getevecfv (vkl(1, ik), vgkl(:, :, :, ik), evecfv)
         Call getevecsv (vkl(1, ik), evecsv)
! find the matching coefficients
         Do ispn = 1, nspnfv
            Call match(ngk(ispn,ik), gkc(:,ispn,ik), &
           & tpgkc(:,:,ispn,ik), sfacgk(:,:,ispn,ik), apwalm(:,:,:,:,ispn))
         End Do
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
! generate the density matrix
               Call gendmat (.False., .False., 0, lmax, is, ia, &
              & ngk(:, ik), apwalm, evecfv, evecsv, lmmax, dmat)
! convert (l,m) part to an irreducible representation if required
               Do ist = 1, nstsv
                  Do ispn = 1, nspinor
                     Do jspn = 1, nspinor
                        Call zgemm ('N', 'N', lmmax, lmmax, lmmax, &
                       & zone, ulm(:, :, ias), lmmax, dmat(:,:,ispn,jspn,ist), &
                       & lmmax, zzero, a, lmmax)
                        Call zgemm ('N', 'C', lmmax, lmmax, lmmax, &
                       & zone, a, lmmax, ulm(:, :, ias), lmmax, &
                       & zzero, dmat(:, :, ispn, jspn, ist), lmmax)
                     End Do
                  End Do
               End Do
! spin rotate the density matrices to desired spin-quantisation axis
               If (associated(input%groundstate%spin) .And. ( .Not. tsqaz)) Then
                  Do ist = 1, nstsv
                     Do lm = 1, lmmax
                        dm1 (:, :) = dmat (lm, lm, :, :, ist)
                        Call z2mm (su2, dm1, dm2)
                        Call z2mmct (dm2, su2, dm1)
                        dmat (lm, lm, :, :, ist) = dm1 (:, :)
                     End Do
                  End Do
               End If
! determine the band characters from the density matrix
               Do ist = 1, nstsv
                  Do ispn = 1, nspinor
                     Do lm = 1, lmmax
                        t1 = dble (dmat(lm, lm, ispn, ispn, ist))
                        bc (lm, ispn, ias, ist, ik) = real (t1)
                     End Do
                  End Do
               End Do
            End Do
         End Do
! compute the spin density matrices of the second-variational states
         Call gensdmat (evecsv, sdmat(:, :, :, ik))
! spin rotate the density matrices to desired spin-quantisation axis
         If (associated(input%groundstate%spin) .And. ( .Not. tsqaz)) &
        & Then
            Do ist = 1, nstsv
               Call z2mm (su2, sdmat(:, :, ist, ik), dm1)
               Call z2mmct (dm1, su2, sdmat(:, :, ist, ik))
            End Do
         End If

         End If

      End Do

! generate energy grid
      dw = (input%properties%dos%winddos(2)-input%properties%dos%winddos(1)) / dble (input%properties%dos%nwdos - 1)
      Do iw = 1, input%properties%dos%nwdos
         w (iw) = dw * dble (iw-1) + input%properties%dos%winddos (1)
      End Do
! number of subdivisions used for interpolation
      nsk (:) = Max(input%properties%dos%ngrdos/input%groundstate%ngridk(:), 1)

!-------------------------------------!
!     calculate and output total DOS  !
!-------------------------------------!

if (rank==0) then

      Open (50, File='TDOS.OUT', Action='WRITE', Form='FORMATTED')
      Call xml_OpenFile ("dos.xml", xf, replace=.True., pretty_print=.True.)
      Call xml_NewElement (xf, "dos")
      Call xml_NewElement (xf, "title")
      Call xml_AddCharacters (xf, trim(input%title))
      Call xml_endElement (xf, "title")
      Call xml_NewElement (xf, "axis")
      Call xml_AddAttribute (xf, "label", 'Energy')
      Call xml_AddAttribute (xf, "unit", 'Hartree')
      Call xml_endElement (xf, "axis")
      Call xml_NewElement (xf, "axis")
      Call xml_AddAttribute (xf, "label", 'DOS')
      Call xml_AddAttribute (xf, "unit", 'states/Hartree/unit cell')
      Call xml_endElement (xf, "axis")
      Call xml_NewElement (xf, "totaldos")
      Do ispn = 1, nspinor
         Call xml_NewElement (xf, "diagram")
         Call xml_AddAttribute (xf, "type", "totaldos")
         Write (buffer,*) ispn
         Call xml_AddAttribute (xf, "nspin", trim(adjustl(buffer)))
         If (ispn .Eq. 1) Then
            t1 = 1.d0
         Else
            t1 = - 1.d0
         End If
         Do ik = 1, nkpt
            Do ist = 1, nstsv
! subtract the Fermi energy
               e (ist, ik, ispn) = evalsv (ist, ik) - efermi
! correction for scissors operator
               If (e(ist, ik, ispn) .Gt. 0.d0) e (ist, ik, ispn) = &
              & e (ist, ik, ispn) + input%properties%dos%scissor
! use diagonal of spin density matrix for weight
               f (ist, ik) = 1.d0
            End Do
         End Do
! BZ integration
         if( input%properties%dos%newint) then
           Call brzint_new (input%properties%dos%nsmdos, &
          & input%groundstate%ngridk, nsk, ikmap, &
          & input%properties%dos%nwdos, input%properties%dos%winddos, nstsv, nstsv, &
          & e(:,:,ispn), f, g(:,ispn))
         else
           Call brzint (input%properties%dos%nsmdos, &
          & input%groundstate%ngridk, nsk, ikmap, &
          & input%properties%dos%nwdos, input%properties%dos%winddos, nstsv, nstsv, &
          & e(:,:,ispn), f, g(:,ispn))
         end if
! multiply by the maximum occupancy (spin-polarised: 1, unpolarised: 2)
         g (:, ispn) = occmax * g (:, ispn)
         Do iw = 1, input%properties%dos%nwdos
            Write (50, '(2G18.10)') w (iw), t1 * g (iw, ispn)
            Call xml_NewElement (xf, "point")
            Write (buffer, '(G18.10)') w (iw)
            Call xml_AddAttribute (xf, "e", trim(adjustl(buffer)))
            Write (buffer, '(G18.10)') g (iw, ispn)
            Call xml_AddAttribute (xf, "dos", trim(adjustl(buffer)))
            Call xml_endElement (xf, "point")
         End Do
         Write (50, '("     ")')
         Call xml_endElement (xf, "diagram")
      End Do
      Close (50)
      Call xml_endElement (xf, "totaldos")
      
!-------------------!
!     Joint DOS     !
!-------------------!
      If (input%properties%dos%jdos) Then
        open(50, File='TJDOS.OUT', Action='WRITE', Form='FORMATTED')
        open(51, File='JDOS.OUT', Action='WRITE', Form='FORMATTED')
        allocate(edif(nstsv,nstsv,nkpt))
        do ispn = 1, nspinor
          if (ispn .Eq. 1) then
            t1 = 1.d0
          else
            t1 = -1.d0
          end if
          edif(:,:,:) = 0.d0
          do ik = 1, nkpt
            write(*,'(i5,3f13.6)') ik, vkl(:,ik)
            n = 0
            do ist = nstsv, 1, -1
            if (evalsv(ist,ik)<=efermi) then
              n = n+1
              m = 0
              do jst = 1, nstsv
              if (evalsv(jst,ik)>efermi) then
                m = m+1
                edif(n,m,ik) = evalsv(jst,ik)-evalsv(ist,ik)+input%properties%dos%scissor
              end if
              end do
              write(*,'(2i5,100f23.16)') ist, n, edif( n, 1:5, ik)
            end if
            end do 
          end do ! ik
          write(*,*) m, n
          ! State dependent JDOS 
          allocate(ej(m,nkpt))
          allocate(fj(m,nkpt))
          fj(:,:) = 1.d0
          do ist = 1, n
            ej(:,:) = edif(ist,1:m,:)
            if( input%properties%dos%newint) then
              call brzint_new(input%properties%dos%nsmdos, &
              &           input%groundstate%ngridk, nsk, ikmap, &
              &           input%properties%dos%nwdos, input%properties%dos%winddos, &
              &           m, m, ej, fj, gp)
            else
              call brzint(input%properties%dos%nsmdos, &
              &           input%groundstate%ngridk, nsk, ikmap, &
              &           input%properties%dos%nwdos, input%properties%dos%winddos, &
              &           m, m, ej, fj, gp)
            end if
            gp(:) = occmax*gp(:)
            do iw = 1, input%properties%dos%nwdos
              if (dabs(w(iw))>1.d-4) then
                write(51,'(3G18.10)') w(iw), t1*gp(iw)/(w(iw)*w(iw))/dble(n*m), t1*gp(iw)
              else
                write(51,'(3G18.10)') w(iw), 0.d0, t1*gp(iw)
              end if  
            end do
            write(51,'("     ")')
          end do ! ist
          deallocate(ej,fj)
          
! total JDOS (sum over all transitions)
          allocate(ej(n*m,nkpt))
          allocate(fj(n*m,nkpt))
          fj(:,:) = 1.d0
          i = 0
          do ist = 1, n
          do jst = 1, m
            i = i+1
            ej(i,:) = edif(ist,jst,:)
          end do
          end do
          write(*,*) input%properties%dos%nsmdos, nsk, input%properties%dos%nwdos, input%properties%dos%winddos
          if( input%properties%dos%newint) then
            call brzint_new(input%properties%dos%nsmdos, &
            &           input%groundstate%ngridk, nsk, ikmap, &
            &           input%properties%dos%nwdos, input%properties%dos%winddos, &
            &           n*m, n*m, ej, fj, gp)
          else
            call brzint(input%properties%dos%nsmdos, &
            &           input%groundstate%ngridk, nsk, ikmap, &
            &           input%properties%dos%nwdos, input%properties%dos%winddos, &
            &           n*m, n*m, ej, fj, gp)
          end if
          gp(:) = occmax*gp(:)
          do iw = 1, input%properties%dos%nwdos
            if (dabs(w(iw))>1.d-4) then
              write(50,'(3G18.10)') w(iw), t1*gp(iw)/(w(iw)*w(iw))/dble(n*m), t1*gp(iw)
            else
              write(50,'(3G18.10)') w(iw), 0.d0, t1*gp(iw)
            end if  
          end do
          write(50, '("     ")')
          deallocate(ej,fj)
          
        end do ! ispn
        deallocate(edif)
        close(50)
        close(51)
      end if ! JDOS

!----------------------------!
!     output partial DOS     !
!----------------------------!

      If (input%properties%dos%lmirep) Then
!
        Do is = 1, nspecies
          Do ia = 1, natoms (is)
            Call xml_NewElement (xf, "partialdos")
            Call xml_AddAttribute (xf, "type", "partial")
            Call xml_AddAttribute (xf, "speciessym", &
           &  trim(adjustl(input%structure%speciesarray(is)%species%chemicalSymbol)))
            Write (buffer,*) is
            Call xml_AddAttribute (xf, "speciesrn", trim(adjustl(buffer)))
            Write (buffer,*) ia
            Call xml_AddAttribute (xf, "atom", trim(adjustl(buffer)))
            ias = idxas (ia, is)
            Write (fname, '("PDOS_S", I2.2, "_A", I4.4, ".OUT")') is, ia
            Open (50, File=trim(fname), Action='WRITE', Form='FORMATTED')
            Do ispn = 1, nspinor
               If (ispn .Eq. 1) Then
                  t1 = 1.d0
               Else
                  t1 = - 1.d0
               End If
               Do l = 0, lmax
                  Do m = - l, l
                     lm = idxlm (l, m)
                     Do ik = 1, nkpt
                        Do ist = 1, nstsv
                           f (ist, ik) = bc (lm, ispn, ias, ist, ik)
                        End Do
                     End Do
                     if( input%properties%dos%newint) then
                       Call brzint_new (input%properties%dos%nsmdos, &
                      &  input%groundstate%ngridk, nsk, ikmap, &
                      &  input%properties%dos%nwdos, input%properties%dos%winddos, &
                      &  nstsv, nstsv, e(:, :, ispn), f, gp)
                     else
                       Call brzint (input%properties%dos%nsmdos, &
                      &  input%groundstate%ngridk, nsk, ikmap, &
                      &  input%properties%dos%nwdos, input%properties%dos%winddos, &
                      &  nstsv, nstsv, e(:, :, ispn), f, gp)
                     end if
                     gp (:) = occmax * gp (:)
                     Call xml_NewElement (xf, "diagram")
                     Write (buffer,*) ispn
                     Call xml_AddAttribute (xf, "nspin", trim(adjustl(buffer)))
                     Write (buffer,*) l
            	     Call xml_AddAttribute (xf, "l", trim(adjustl(buffer)))
                     Write (buffer,*) m
            	     Call xml_AddAttribute (xf, "m", trim(adjustl(buffer)))
                     Do iw = 1, input%properties%dos%nwdos
                        Call xml_NewElement (xf, "point")
                        Write (buffer, '(G18.10)') w (iw)
                        Call xml_AddAttribute (xf, "e", trim(adjustl(buffer)))
                        Write (buffer, '(G18.10)') gp (iw)
                        Call xml_AddAttribute (xf, "dos", trim(adjustl(buffer)))
                        Call xml_endElement (xf, "point")
                        Write (50, '(2G18.10)') w (iw), t1 * gp (iw)
                        g (iw, ispn) = g (iw, ispn) - gp (iw)
                     End Do
                     Write (50, '("     ")')
                   Call xml_endElement (xf, "diagram")
                  End Do
               End Do
            End Do
            Close (50)
            Call xml_endElement (xf, "partialdos")
          End Do
        End Do
!
        Open (50, File='ELMIREP.OUT', Action='WRITE', Form='FORMATTED')
        Call xml_NewElement (xf, "limrep")
        Do is = 1, nspecies
          Call xml_NewElement (xf, "species")
          Call xml_AddAttribute (xf, "speciessym", trim(adjustl(input%structure%speciesarray(is)%species%chemicalSymbol)))
          Do ia = 1, natoms (is)
            Call xml_NewElement (xf, "atom")
            ias = idxas (ia, is)
            Write (50,*)
            Write (50, '("Species : ", I4, " (", A, "), atom : ", I4)') &
           &  is, trim (spsymb(is)), ia
            Do l = 0, lmax
               Do m = - l, l
                  lm = idxlm (l, m)
                  Call xml_NewElement (xf, "orb")
                  Write (50, '(" l = ", I2, ", m = ", I2, ", lm= ", I3, &
                 & " : ", G18.10)') l, m, lm, elm (lm, ias)
                  Write (buffer,*) l
                  Call xml_AddAttribute (xf, "l", trim(adjustl(buffer)))
                  Write (buffer,*) m
                  Call xml_AddAttribute (xf, "m", trim(adjustl(buffer)))
                  Write (buffer,*) lm
                  Call xml_AddAttribute (xf, "lm", trim(adjustl(buffer)))
                  Write (buffer, '(G18.10)') elm (lm, ias)
                  Call xml_AddAttribute (xf, "elm", trim(adjustl(buffer)))
                  Call xml_endElement (xf, "orb")
               End Do
            End Do
            Call xml_endElement (xf, "atom")
          End Do
          Call xml_endElement (xf, "species")
        End Do
        Call xml_endElement (xf, "limrep")
        Close (50)

!---------------------------------!
!     output interstitial DOS     !
!---------------------------------!

        Call xml_NewElement (xf, "interstitialdos")
        Open (50, File='IDOS.OUT', Action='WRITE', Form='FORMATTED')
        Do ispn = 1, nspinor
          Call xml_NewElement (xf, "diagram")
          Call xml_AddAttribute (xf, "type", "interstitial")
          Write (buffer,*) ispn
          Call xml_AddAttribute (xf, "nspin", trim(adjustl(buffer)))
          If (ispn .Eq. 1) Then
            t1 = 1.d0
          Else
            t1 = - 1.d0
          End If
          Do iw = 1, input%properties%dos%nwdos
            Call xml_NewElement (xf, "point")
            Write (buffer, '(G18.10)') w (iw)
            Call xml_AddAttribute (xf, "e", trim(adjustl(buffer)))
            Write (buffer, '(G18.10)') g (iw, ispn)
            Call xml_AddAttribute (xf, "dos", trim(adjustl(buffer)))
            Call xml_endElement (xf, "point")
            Write (50, '(2G18.10)') w (iw), t1 * g (iw, ispn)
          End Do
          Call xml_endElement (xf, "diagram")
        End Do
        Call xml_endElement (xf, "interstitialdos")

      End If ! lmirep

! close files
      Close (50)
      Call xml_endElement (xf, "dos")
      Call xml_close (xf)
      
      Deallocate (e, f, w, g, gp, bc)
      If (input%properties%dos%lmirep) deallocate (elm, ulm, a)
      Deallocate (dmat, sdmat, apwalm, evecfv, evecsv)

!-----------------------------------
! Screen info
!-----------------------------------
      Write(*,*)
      Write(*, '("Info(dos):")')
      Write(*,*)
      Write(*, '("   Total density of states written to TDOS.OUT")')
      Write(*,*)
      If (input%properties%dos%lmirep) Then
         Write(*, '("   Partial density of states written to PDOS_Sss_Aaaaa.OUT")')
         Write(*, '("   for all species and atoms")')
         Write(*,*)
         Write(*, '("   Eigenvalues of a random matrix in the (l, m) basis symmetrised")')
         Write(*, '("   with the site symmetries written to ELMIREP.OUT for all")')
         Write(*, '("   species and atoms. Degenerate eigenvalues correspond to")')
         Write(*, '("   irreducible representations of each site symmetry group")')
         Write(*,*)
         Write(*, '("   Interstitial density of states written to IDOS.OUT")')
         Write(*,*)
      End If
      If (input%properties%dos%jdos) Then
         Write(*, '("   Joint density of states written to JDOS.OUT")')
         Write(*,*)
      End If
      Write(*, '("   Fermi energy is at zero in plot")')
      Write(*,*)
      Write(*, '("   DOS units are states/Hartree/unit cell")')
      Write(*,*)

end if ! rank

      Return
End Subroutine
!EOC
