!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: rhoplot
! !INTERFACE:
!
Subroutine rhoplot
! !USES:
      Use modinput
      Use modmain
      use modplotlabels
      use modmpi, only : rank

! !DESCRIPTION:
!   Outputs the charge density and the charge density gradients (modulus)
!   read in from {\tt STATE.OUT}, for 1D, 2D or 3D
!   plotting.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!   Density gradients are added, March 2011 (DIN)
!EOP
!BOC
      Implicit None
! Local variables for gradient calculation
      real(8), allocatable :: rhomtsp(:,:,:,:)
      real(8), allocatable :: rhoirsp(:,:)
      Integer :: i, j, k
! initialise universal variables
      type(plotlabels), pointer :: labels

!---------------------------------------------
! Read density (and magnetization) from file
!---------------------------------------------
      Call init0
! read density from file
      Call readstate
! visualize only valence density      
      if (input%properties%chargedensityplot%nocore) Then
        call gencore
        call removerhocr
      end if
     
!----------------------------------------
! Calculate spin-up and -down densities
!----------------------------------------
      if (associated(input%groundstate%spin)) then
        if (allocated(rhomtsp)) deallocate (rhomtsp)
        allocate(rhomtsp(lmmaxvr,nrmtmax,natmtot,2))
        if (allocated(rhoirsp)) deallocate (rhoirsp)
        allocate(rhoirsp(ngrtot,2))
        call calcrhospinpol(rhomtsp,rhoirsp)
      end if

!---------------------------------      
! write the density plot to file
!---------------------------------
      ! 1D plot
      if (associated(input%properties%chargedensityplot%plot1d)) Then
        if (associated(input%groundstate%spin)) then
if (.false.) then
          labels=>create_plotlablels("RHOUP1D","RHOUP1D",1)
          call set_plotlabel_axis(labels,1,"Distance","a_0","graceunit")
          call set_plotlabel_axis(labels,2,"Density","???","graceunit")
          call plot1d(labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
          &           rhomtsp(:,:,:,1), rhoirsp(:,1), &
          &           input%properties%chargedensityplot%plot1d)
          call destroy_plotlablels(labels)
          if (rank==0) then
            write(*,*)
            write(*, '("Info(rhoplot):")')
            write(*, '(" 1D spin-up density plot written to RHOUP1D.xml")')
          end if
          labels=>create_plotlablels("RHODN1D","RHODN1D",1)
          call set_plotlabel_axis(labels,1,"Distance","a_0","graceunit")
          call set_plotlabel_axis(labels,2,"Density","???","graceunit")
          call plot1d(labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
          &           rhomtsp(:,:,:,2), rhoirsp(:,2), &
          &           input%properties%chargedensityplot%plot1d)
          call destroy_plotlablels(labels)
          if (rank==0) then
            write(*, '(" 1D spin-down density plot written to RHODN1D.xml")')
          end if
end if
          labels=>create_plotlablels("RHO1D","RHO1D",1)
          call set_plotlabel_axis(labels,1,"Distance","a_0","graceunit")
          call set_plotlabel_axis(labels,2,"Density","???","graceunit")
          call plot1d(labels, 2, input%groundstate%lmaxvr, lmmaxvr, &
          &           rhomtsp, rhoirsp, &
          &           input%properties%chargedensityplot%plot1d)
          call destroy_plotlablels(labels)
          if (rank==0) then
            write(*, '("Info(rhoplot):")')
            write(*, '(" 1D spin-dependent density plot written to RHO1D.xml")')          
          end if
        else
          labels=>create_plotlablels("RHO1D","RHO1D",1)
          call set_plotlabel_axis(labels,1,"Distance","a_0","graceunit")
          call set_plotlabel_axis(labels,2,"Density","???","graceunit")
          call plot1d (labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
          &            rhomt, rhoir, input%properties%chargedensityplot%plot1d)
          call destroy_plotlablels(labels)
          if (rank==0) then
            write(*, '("Info(rhoplot):")')
            write(*, '(" 1D density plot written to RHO1D.xml")')
          end if
        end if
      end if
      
      ! 2D plot
      If (associated(input%properties%chargedensityplot%plot2d)) Then
        if (associated(input%groundstate%spin)) then
if (.false.) then        
          labels=>create_plotlablels("RHOUP2D","RHOUP2D",2)
	        call set_plotlabel_axis(labels,1,"a","1","graceunit")
	        call set_plotlabel_axis(labels,2,"b","1","graceunit")
	        call set_plotlabel_axis(labels,3,"Density","???","graceunit")
          call plot2d(labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
          &           rhomtsp(:,:,:,1), rhoirsp(:,1), &
          &           input%properties%chargedensityplot%plot2d)
          call destroy_plotlablels(labels)
          if (rank==0) then
            write(*, '("Info(rhoplot):")')
            Write(*, '(" 2D spin-up density plot written to RHOUP2D.xml")')
          end if
	        labels=>create_plotlablels("RHODN2D","RHODN2D",2)
	        call set_plotlabel_axis(labels,1,"a","1","graceunit")
	        call set_plotlabel_axis(labels,2,"b","1","graceunit")
	        call set_plotlabel_axis(labels,3,"Density","???","graceunit")
          call plot2d(labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
          &           rhomtsp(:,:,:,2), rhoirsp(:,2), &
          &           input%properties%chargedensityplot%plot2d)
          call destroy_plotlablels(labels)
          if (rank==0) then
            Write(*, '(" 2D spin-down density plot written to RHODN2D.xml")')
          end if
end if
	        labels=>create_plotlablels("RHO2D","RHO2D",2)
	        call set_plotlabel_axis(labels,1,"a","1","graceunit")
	        call set_plotlabel_axis(labels,2,"b","1","graceunit")
	        call set_plotlabel_axis(labels,3,"Density","???","graceunit")
          call plot2d(labels, 2, input%groundstate%lmaxvr, lmmaxvr, &
          &           rhomtsp, rhoirsp, &
          &           input%properties%chargedensityplot%plot2d)
          call destroy_plotlablels(labels)
          if (rank==0) then
            write(*, '("Info(rhoplot):")')
            Write(*, '(" 2D spin-dependent density plot written to RHO2D.xml")')
          end if
        else
          labels=>create_plotlablels("RHO2D","RHO2D",2)
	        call set_plotlabel_axis(labels,1,"a","1","graceunit")
	        call set_plotlabel_axis(labels,2,"b","1","graceunit")
	        call set_plotlabel_axis(labels,3,"Density","???","graceunit")
          call plot2d(labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
          &           rhomt, rhoir, input%properties%chargedensityplot%plot2d)
          call destroy_plotlablels(labels)
          if (rank==0) then
            write(*, '("Info(rhoplot):")')
            write(*, '(" 2D density plot written to RHO2D.xml")')
          end if
        end if
      End If
      
      ! 3D plot
      If (associated(input%properties%chargedensityplot%plot3d)) Then
        if (associated(input%groundstate%spin)) then
if (.false.) then        
          labels=>create_plotlablels("RHOUP3D","RHOUP3D",3)
          call set_plotlabel_axis(labels,1,"a","1","graceunit")
          call set_plotlabel_axis(labels,2,"b","1","graceunit")
          call set_plotlabel_axis(labels,3,"b","1","graceunit")
          call set_plotlabel_axis(labels,4,"Density","???","graceunit")
          Call plot3d(labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
          &           rhomtsp(:,:,:,1), rhoirsp(:,1), &
          &           input%properties%chargedensityplot%plot3d)
          call destroy_plotlablels(labels)
          if (rank==0) then
            write(*, '("Info(rhoplot):")')
            Write(*, '(" 3D spin-up density plot written to RHOUP3D.xml")')
          end if
          labels=>create_plotlablels("RHODN3D","RHODN3D",3)
          call set_plotlabel_axis(labels,1,"a","1","graceunit")
          call set_plotlabel_axis(labels,2,"b","1","graceunit")
          call set_plotlabel_axis(labels,3,"b","1","graceunit")
          call set_plotlabel_axis(labels,4,"Density","???","graceunit")
          Call plot3d(labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
          &           rhomtsp(:,:,:,2), rhoirsp(:,2), &
          &           input%properties%chargedensityplot%plot3d)
          call destroy_plotlablels(labels)
          if (rank==0) then
            Write(*, '(" 3D spin-down density plot written to RHODN3D.xml")')
          end if
end if
          labels=>create_plotlablels("RHO3D","RHO3D",3)
          call set_plotlabel_axis(labels,1,"a","1","graceunit")
          call set_plotlabel_axis(labels,2,"b","1","graceunit")
          call set_plotlabel_axis(labels,3,"b","1","graceunit")
          call set_plotlabel_axis(labels,4,"Density","???","graceunit")
          Call plot3d(labels, 2, input%groundstate%lmaxvr, lmmaxvr, &
          &           rhomtsp, rhoirsp, &
          &           input%properties%chargedensityplot%plot3d)
          call destroy_plotlablels(labels)
          if (rank==0) then
            write(*, '("Info(rhoplot):")')
            write(*, '(" 3D spin-dependent density plot written to RHO3D.xml")')
          end if
        else
          labels=>create_plotlablels("RHO3D","RHO3D",3)
          call set_plotlabel_axis(labels,1,"a","1","graceunit")
          call set_plotlabel_axis(labels,2,"b","1","graceunit")
          call set_plotlabel_axis(labels,3,"b","1","graceunit")
          call set_plotlabel_axis(labels,4,"Density","???","graceunit")
          Call plot3d(labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
          &           rhomt, rhoir, input%properties%chargedensityplot%plot3d)
          !Call plot3d_ir(labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
          !&              rhomt, rhoir, input%properties%chargedensityplot%plot3d)
          call destroy_plotlablels(labels)
          if (rank==0) then
            write(*, '("Info(rhoplot):")')
            write(*, '(" 3D density plot written to RHO3D.xml")')
          end if
        end if
      End If
      if (rank==0) write(*,*)
      Return
      
contains

    subroutine calcrhospinpol(rhomtsp,rhoirsp)
        use modinput
        use modmain
        implicit none
        ! output
        real(8), intent(Out) :: rhomtsp(lmmaxvr,nrmtmax,natmtot,2)
        real(8), intent(Out) :: rhoirsp(ngrtot,2)
        ! local variables
        integer :: i, n, is, ia, ias, ir, nr, idm, lm
        real(8) :: t1, t2, bext(3)
        real(8), allocatable :: rho(:), rhoup(:), rhodn(:), mag(:,:)

        !------------
        ! muffin-tin 
        !------------
        rhomtsp(:,:,:,:) = 0.d0
        n = lmmaxvr*nrmtmax
        allocate(rho(n),rhoup(n),rhodn(n),mag(n,3))
        do is = 1, nspecies
          nr = nrmt(is)
          n = lmmaxvr*nr
          do ia = 1, natoms(is)
            ias = idxas(ia,is)
            ! compute the density in spherical coordinates
            call dgemm('N','N',lmmaxvr,nr,lmmaxvr, &
            &          1.d0,rbshtvr,lmmaxvr, &
            &          rhomt(:,:,ias),lmmaxvr, &
            &          0.d0,rho,lmmaxvr)
            ! magnetisation in spherical coordinates
            do idm = 1, ndmag
              call dgemm('N','N',lmmaxvr,nr,lmmaxvr, &
              &          1.d0,rbshtvr,lmmaxvr, &
              &          magmt(:,:,ias,idm),lmmaxvr, &
              &          0.d0,mag(:,idm),lmmaxvr)
            end do
            if (ncmag) then
              ! non-collinear (use Kubler's trick)
              bext(:) = input%groundstate%spin%bfieldc(:)+ &
              &         input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:)
              do i = 1, n
                ! compute rhoup=(rho+sgn(m.B_ext)|m|)/2 and rhodn=(rho-sgn(m.B_ext)|m|)/2
                t1 = sqrt(mag(i,1)**2+mag(i,2)**2+mag(i,3)**2)
                if (xcgrad.ne.0) then
                  t2 = mag(i,1)*bext(1)+mag(i,2)*bext(2)+mag(i,3)*bext(3)
                  if (t2<0.d0) t1=-t1
                end if
                rhoup(i) = 0.5d0*(rho(i)+t1)
                rhodn(i) = 0.5d0*(rho(i)-t1)
              end do
            else
              ! collinear
              do i = 1, n
                ! compute rhoup=(rho+m_z)/2 and rhodn=(rho-m_z)/2
                rhoup(i) = 0.5d0*(rho(i)+mag(i,1))
                rhodn(i) = 0.5d0*(rho(i)-mag(i,1))
              end do
            end if ! ncmag
            ! convert density into spherical harmonics
            call dgemm('N','N',lmmaxvr,nr,lmmaxvr, &
            &          1.d0,rfshtvr,lmmaxvr, &
            &          rhoup,lmmaxvr, &
            &          0.d0,rhomtsp(:,:,ias,1),lmmaxvr)
            call dgemm('N','N',lmmaxvr,nr,lmmaxvr, &
            &          1.d0,rfshtvr,lmmaxvr, &
            &          rhodn,lmmaxvr, &
            &          0.d0,rhomtsp(:,:,ias,2),lmmaxvr)
          end do ! ia
        end do ! is
        deallocate(rho,rhoup,rhodn,mag)
       
        !--------------
        ! interstitial 
        !--------------
        rhoirsp(:,:) = 0.d0
        if (ncmag) then
          ! non-collinear
          do ir = 1, ngrtot
            t1 = sqrt(magir(ir,1)**2+magir(ir,2)**2+magir(ir,3)**2)
            if (xcgrad.ne.0) then
              t2 = magir(ir,1)*input%groundstate%spin%bfieldc(1)+ &
              &    magir(ir,2)*input%groundstate%spin%bfieldc(2)+ &
              &    magir(ir,3)*input%groundstate%spin%bfieldc(3)
              if (t2<0.d0) t1=-t1
            end if
            rhoirsp(ir,1) = 0.5d0*(rhoir(ir)+t1)
            rhoirsp(ir,2) = 0.5d0*(rhoir(ir)-t1)
          end do
        else
          ! collinear
          do ir = 1, ngrtot
            rhoirsp(ir,1) = 0.5d0*(rhoir(ir)+magir(ir,1))
            rhoirsp(ir,2) = 0.5d0*(rhoir(ir)-magir(ir,1))
          end do
        end if ! ncmag
        return
    end subroutine
      
End Subroutine
!EOC
