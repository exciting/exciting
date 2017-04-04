!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOI
! !TITLE: The XS/EXCITING Code (eXited States) Manual \\ Version 0.9
! !AUTHORS: S. Sagmeister and C. Ambrosch-Draxl
! !AFFILIATION:
! !INTRODUCTION: Introduction
!   Welcome to the {\sf XS/EXCITING} code developers' manual.
!   This manual is supposed to collect the routines and modules belonging
!   exclusively to the excited states (TDDFT and BSE) part into one document.
!   The content of this manual is partially taken from the author's PhD thesis.
!   \\\&
!   S. Sagmeister\&
!   Leoben, August 2008
!
!EOI
!
Module modxsmain
      Use mod_misc
      Implicit None
      Integer :: nxstasks
      Character (256) :: xstasks (maxtasks)
End Module
!
!
Subroutine xsmain (plan)
      Use modmain
      Use modinput
      Use modmpi
#ifdef TETRA      
      Use modtetra
#endif
      Use modxs
      use mod_exciton_wf
      Implicit None
      Type (plan_type) :: plan
      Integer :: i
  ! initialization
!
  ! task selection
      Do i = 1, size (plan%doonlyarray)
         task = plan%doonlyarray(i)%doonly%tasknumber
         Call xsinit
         Select Case (task)
         Case (23)
     ! estimate bandgap from regular grid
            Call writebandgapgrid
         Case (700)
     ! estimate disk-space, cpu-time and memory
            Call xsestimate
         Case (701)
     ! test timing
            Call xstiming
         Case (301)
     ! generate eigenvectors, eigenvalues, occupancies and MT-coefficients
     ! for q-point set
            Call xsgeneigvec
         Case (310)
#ifdef TETRA          
     ! calculate weights for tetrahedron method
            Call tetcalccw
#else
            ! added by DIN
            write(*,*) 'Tetrahedron method for XS is disabled!'
            write(*,*) 'Check -DTETRA option in make.inc' 
            stop
#endif            
         Case (320)
     ! parallel version of momentum matrix elements
			 if (input%xs%BSE%xas) then
				Call xasinit
				Call writepmatxs
				Call xasfinit
			else
				Call writepmatxs
			end if
         Case (321)
     ! ASCII output of momentum matrix elements
        if (input%xs%BSE%xas) then
          Call xasinit
          Call writepmatasc
          Call xasfinit
        else
          Call writepmatasc
        end if
         Case (322)
     ! convert momentum matrix elements file to old format
            Call pmatxs2orig
         Case (330)
     ! calculate matrix elements of exponential expression (band combs)
            Call writeemat
         Case (331)
     ! ASCII output of matrix elements of exponential expression
            Call writeematasc
         Case (335)
     ! calculate matrix elements of the plane wave (simple version for checking)
            Call writepwmat
         Case (339)
     ! check relation between matr. el. of exp. and mom. matr. el.
            Call emattest
         Case (340)
     ! Kohn Sham response function
            Call df
         Case (341)
     ! ASCII output of Kohn Sham response function
            Call x0toasc
         Case (342)
     ! binary output of Kohn Sham response function
            Call x0tobin
         Case (350)
     ! inverse of dielectric function - solve Dyson equation for xc-kernel
            Call idf
         Case (398)
     ! check ALDA kernel
            Call fxc_alda_check
         Case (401)
     ! generate eigenvectors, eigenvalues, occupancies and APW MT coefficients
     ! for screening and BSE(-kernel)
            Call scrgeneigvec
         Case (410)
#ifdef TETRA         
     ! calculate weights for tetrahedron method (screening)
            Call scrtetcalccw
#else
            ! added by DIN
            write(*,*) 'Tetrahedron method for XS is disabled!'
            write(*,*) 'Check -DTETRA option in make.inc' 
            stop
#endif            
         Case (420)
     ! momentum matrix elements for screening
            Call scrwritepmat
         Case (430)
     ! RPA screening
            Call screen
         Case (440)
     ! screened Coulomb interaction
			if (input%xs%BSE%xas) then
				Call xasinit
				Call xas_scrcoulint
				Call xasfinit
			else
				Call scrcoulint
			end if
         Case (441)
     ! exchange Coulomb interaction
			 if (input%xs%BSE%xas) then
				Call xasinit
				Call xas_exccoulint
				Call xasfinit
			else
				Call exccoulint
			end if
         Case (445)
     ! Bethe-Salpeter equation
             if (input%xs%BSE%xas) then
				Call xasinit
				Call xas
				Call xasfinit
			else
	            Call BSE
	        end if
         Case (446)
     ! regenerate BSE spectrum from exciton output

            if (input%xs%BSE%xas) then
				call xasinit
				call bsegenspec
				call xasfinit
			else
				Call bsegenspec
			end if

            Call bsegenspec
         Case (447)
     ! ASCII output of BSE eigenvectors
            Call writebevec
         Case (450)
     ! BSE-kernel
            Call kernxc_bse
         Case (451)
     ! BSE-kernel, simple version, very slow
            Call kernxc_bse3
         Case (499)
     ! call to test xs-routine
            Call testxs
         case (999)
            call testmain
         case (710)
			if (input%xs%BSE%xas) then
				call xasinit
				call plot_excitonWavefunction
				call xasfinit
			else
				call plot_excitonWavefunction
			end if
         Case Default
            Write (*,*)
            Write (*,*) 'Error(xsmain): task not defined:', task
            Write (*,*)
            Call terminate
         End Select
         Call xsfinit
      End Do
  ! summarize information on run
!
End Subroutine xsmain
