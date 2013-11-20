!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: rhoinit
! !INTERFACE:
!
!
Subroutine rhoinit
! !USES:
      Use modinput
      Use modmain
      use omp_lib
! !DESCRIPTION:
!   Initialises the crystal charge density. Inside the muffin-tins it is set to
!   the spherical atomic density. In the interstitial region it is taken to be
!   constant such that the total charge is correct. Requires that the atomic
!   densities have already been calculated.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
! polynomial order of smooth step function
      Integer, Parameter :: n = 4
      Integer :: lmax, lmmax, l, m, lm, ir, irc
      Integer :: is, ia, ias, ig, ifg
      Real (8) :: x, t1, t2, jlgr01
      Complex (8) zt1, zt2, zt3
      Real (8) :: ta,tb, tc,td
! automatic arrays
      Real (8) :: fr (spnrmax), gr (spnrmax), cf (3, spnrmax)
! allocatable arrays
      Real (8), Allocatable :: jlgr (:)
      Real (8), Allocatable :: th (:, :)
      Real (8), Allocatable :: ffacg (:)
      Complex (8), Allocatable :: zfmt (:, :),z2fmt (:)
      Complex (8), Allocatable :: zfft (:)

      integer :: auxgridsize,mtgridsize,lastpoint
      Real (8), Allocatable :: auxgrid(:),auxrho(:),a(:),c(:)
      Real (8) :: rhoder,rhoder2
      real(8), parameter :: threshold=1d-6
      integer, parameter :: PointsPerPeriod=10

      integer:: boundhi,boundlo
      real(8) :: cpu_seconds
      external :: cpu_seconds

! maximum angular momentum for density initialisation
      lmax = 1
      lmmax = (lmax+1) ** 2
! allocate local arrays
!      Allocate (jlgr(0:lmax, nrcmtmax))
      Allocate (jlgr(0:lmax))
      Allocate (ffacg(ngvec))
      Allocate (zfmt(lmmax, nrcmtmax))
      Allocate (z2fmt(lmmax))
      Allocate (zfft(ngrtot))
      allocate(a(nspecies))
      allocate(c(nspecies))

! zero the charge density and magnetisation arrays
      rhomt (:, :, :) = 0.d0
      rhoir (:) = 0.d0
      If (associated(input%groundstate%spin)) Then
         magmt (:, :, :, :) = 0.d0
         magir (:, :) = 0.d0
      End If
      call timesec(ta)
! compute the superposition of all the atomic density tails
      zfft (:) = 0.d0
      Do is = 1, nspecies
         mtgridsize=PointsPerPeriod*int(rmt(is)*input%groundstate%gmaxvr/(2d0*pi))+1
         lastpoint=spnr(is)
         do while ((lastpoint.gt.nrmt(is)).and.(sprho (lastpoint, is).lt.threshold))
           lastpoint=lastpoint-1
         enddo

         auxgridsize=mtgridsize+lastpoint-nrmt(is)+1
         allocate(auxgrid(auxgridsize))
         allocate(auxrho(auxgridsize))
         auxgrid(1)=spr(1,is)
         auxrho(1)=sprho (1, is)

         do ir = 2, mtgridsize
           auxgrid(ir)=rmt(is)*dble(ir-1)/dble(mtgridsize-1)
         enddo
         auxgrid(mtgridsize+1:auxgridsize)=spr(nrmt(is)+1:lastpoint, is)
         Call fderiv (1, spnr(is), spr(:, is), sprho (:, is), gr, cf)
         rhoder=cf(1,nrmt(is))
         rhoder2=cf(2,nrmt(is))*2d0
         a(is)=rhoder/(2d0*rmt(is))
         c(is)=sprho (nrmt(is), is)-a(is)*rmt(is)**2
         do ir = 1, mtgridsize
           auxrho(ir)=a(is)*auxgrid(ir)**2+c(is)
         enddo
         auxrho(mtgridsize+1:auxgridsize)=sprho(nrmt(is)+1:lastpoint, is)
!        do ir = 1, auxgridsize
!          write(*,*) auxgrid(ir),auxrho(ir)
!        enddo

         call timesec(tc)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ig,ir,x,jlgr01,fr,cf,gr)
!$OMP DO
         Do ig = 1, ngvec
           do ir = 1,auxgridsize
              x=gc(ig)*auxgrid(ir)
              Call sbessel (0, x, jlgr01)
              fr (ir) = auxrho(ir) * jlgr01 * auxgrid(ir) ** 2
           enddo
           Call fderiv (-1, auxgridsize, auxgrid, fr, gr, cf)
           ffacg (ig) = (fourpi/omega) * gr (auxgridsize)
         End Do
!$OMP END DO
!$OMP END PARALLEL

         deallocate(auxgrid,auxrho)
         call timesec(td)
         write(*,*) td-tc
         write(*,*)

         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ig = 1, ngvec
               ifg = igfft (ig)
               zfft (ifg) = zfft (ifg) + ffacg (ig) * conjg (sfacg(ig, ias))
            End Do
         End Do
!         call timesec(td)
!         write(*,*) td-tc
!         read(*,*) 
      End Do
      call timesec(tb)
      write(*,*) 'rhoinit, step 1:',tb-ta
      call timesec(ta)
! compute the tails in each muffin-tin

      do ig=1,ngvec
        ffacg(ig)=zfft(igfft (ig))
      enddo
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            zfmt (:, :) = 0.d0
!$OMP PARALLEL DEFAULT(NONE) SHARED(is,ias,rcmt,gc,nrcmt,ngvec,lmax,sfacg,ffacg,zil,ylmg,zfmt) PRIVATE(ig,irc,x,zt1,lm,l,zt2,m,zt3,jlgr,z2fmt)
!$OMP DO
            Do irc = 1, nrcmt (is)
              z2fmt(:)=0d0
              Do ig = 1, ngvec
                x = gc (ig) * rcmt (irc, is)
                if (ig.ne.1) then
                  if (gc(ig)-gc(ig-1).gt.1d-8) Call sbessel (lmax, x, jlgr)
                else
                  Call sbessel (lmax, x, jlgr)
                endif
                zt1 = fourpi * ffacg (ig) * sfacg (ig, ias)
                lm = 0
                Do l = 0, lmax
                  zt2 = zt1 * zil (l)
                  Do m = - l, l
                     lm = lm + 1
                     zt3 = zt2 * conjg (ylmg(lm, ig))
                     z2fmt(lm)=z2fmt(lm)+jlgr(l) * zt3
!                     zfmt (lm, irc) = zfmt (lm, irc) + jlgr (l) * zt3
                  End Do
               End Do
              End Do
              zfmt (:, irc)=z2fmt(:)
            End Do
!$OMP END DO
!$OMP END PARALLEL
            irc = 0
            Do ir = 1, nrmt (is), input%groundstate%lradstep
               irc = irc + 1
               Call ztorflm (lmax, zfmt(:, irc), rhomt(:, ir, ias))
            End Do
         End Do
      End Do

      call timesec(tb)
      write(*,*) 'rhoinit, step 2:',tb-ta
      call timesec(ta)
! convert the density from a coarse to a fine radial mesh
      Call rfmtctof (rhomt)
! add the atomic charge density and the excess charge in each muffin-tin
      t1 = input%groundstate%chgexs / omega
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ir = 1, nrmt (is)
               t2 = (t1+sprho(ir, is)-a(is)*spr(ir,is)**2-c(is)) / y00
               rhomt (1, ir, ias) = rhomt (1, ir, ias) + t2
!               write(*,*) spr(ir,is),rhomt (1, ir, ias)*y00,sprho(ir,is)
            End Do
         End Do
      End Do
!      stop
! interstitial density determined from the atomic tails and excess charge
      Call zfftifc (3, ngrid, 1, zfft)
      Do ir = 1, ngrtot
         rhoir (ir) = dble (zfft(ir)) + t1
      End Do
! compute the total charge
      Call charge
! normalise the density
      Call rhonorm
      Deallocate (jlgr, ffacg, zfmt, zfft,z2fmt,a,c)
      call timesec(tb)
      write(*,*) 'rhoinit, step 3:',tb-ta
!      stop
      Return
End Subroutine
real(8) function cpu_seconds()
implicit none
integer values(8)
call date_and_time(values=values)
cpu_seconds=(values(3)*24+values(5))*3600.d0+values(6)*60.d0+&
  values(7)*1.d0+values(8)/1000.d0
return
end function
!EOC
