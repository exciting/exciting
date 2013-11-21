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
#ifdef USEOMP
      use omp_lib
#endif
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
      Complex (8) zt1, zt2, zt3,yy(4),update(4)
      Real (8) :: ta,tb, tc,td
! automatic arrays
      Real (8) :: fr (spnrmax), gr (spnrmax), cf (3, spnrmax)
! allocatable arrays
      Real (8), Allocatable :: jlgr (:),jj(:,:)
      Real (8), Allocatable :: th (:, :)
      Real (8), Allocatable :: ffacg (:)
      Complex (8), Allocatable :: zfmt (:, :),z2fmt (:,:)
      Complex (8), Allocatable :: zfft (:)

      integer :: auxgridsize,mtgridsize,lastpoint
      Real (8), Allocatable :: auxgrid(:),auxrho(:),a(:),c(:)
      Real (8) :: rhoder,rhoder2, tp(2),r
      real(8), parameter :: threshold=1d-6
      integer, parameter :: PointsPerPeriod=10

      integer:: boundhi,boundlo,middle
      real(8):: jthr,cs,sn,xi
#ifdef USEOMP
      integer:: whichthread,nthreads
#endif
      real(8), parameter :: sixth=1d0/6d0
      real(8), parameter :: third=1d0/3d0
      real(8) :: cpu_seconds
      external :: cpu_seconds

! maximum angular momentum for density initialisation
      lmax = 1
      lmmax = (lmax+1) ** 2
! allocate local arrays
      Allocate (jj(0:lmax, nrcmtmax))
      Allocate (jlgr(0:lmax))
      Allocate (ffacg(ngvec))
      Allocate (zfmt(lmmax, nrcmtmax))
      Allocate (z2fmt(lmmax, nrcmtmax))
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
! The density within the IR region is smoothly continued in the MT region
! by a parabola a*r**2+c. Then, the pseudodensity is smooth and slowly varying
! everywhere, and relatively few grid points are needed for calculating 
! integrals.

! Specifying new grid. It is equally spaced within MT and coincides with the original
! one outside.
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

! Calculating derivatives from splines in order to obtain parameters of a*r**2+c.
         Call fderiv (1, spnr(is), spr(:, is), sprho (:, is), gr, cf)
         rhoder=cf(1,nrmt(is))
         rhoder2=cf(2,nrmt(is))*2d0
         a(is)=rhoder/(2d0*rmt(is))
         c(is)=sprho (nrmt(is), is)-a(is)*rmt(is)**2
! Interior of MT is filled with a*r**2+c
         do ir = 1, mtgridsize
           auxrho(ir)=a(is)*auxgrid(ir)**2+c(is)
         enddo
! Exterior is filled with the actual atomic density.
         auxrho(mtgridsize+1:auxgridsize)=sprho(nrmt(is)+1:lastpoint, is)

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ig,jthr,boundhi,boundlo,ir,x,sn,t1,fr,cf,gr)
!$OMP DO
#endif 
         Do ig = 1, ngvec
! First, we separate x<1d-8 from large x.
! It is done using binary search.
! This block implements the following commented loop.
! See below for the explanations.
!           do ir = 1,auxgridsize
!              x=gc(ig)*auxgrid(ir)
!              Call sbessel (0, x, jlgr01)
!              fr (ir) = auxrho(ir) * jlgr01 * auxgrid(ir) ** 2
!           enddo
               jthr=1d-8/gc (ig)
               boundhi=auxgridsize
               boundlo=1
               if (auxgrid (1).gt.jthr) then
                 boundlo=0
                 boundhi=1
               elseif (auxgrid (auxgridsize).lt.jthr) then
                 boundlo=auxgridsize
                 boundhi=auxgridsize+1
               else
                 do while (boundhi-boundlo.gt.1)
                   ir=(boundhi+boundlo)/2
                   if (auxgrid (ir).gt.jthr) then
                     boundhi=ir
                   else
                     boundlo=ir
                   endif
                 enddo
               endif

               do ir=1,boundlo
                 x=gc (ig)*auxgrid(ir)
                 fr (ir) = auxrho(ir) * (1d0-sixth*x**2) * auxgrid(ir) ** 2
               enddo

               t1=1d0/gc (ig)
               do ir=boundhi,auxgridsize
                 x=gc (ig)*auxgrid(ir)
                 sn=sin(x)
                 fr (ir) = auxrho(ir) * sn*t1 * auxgrid(ir) 
               enddo
! End of block
! Uncomment the following block for the trapezoid rule. 
!           gr(auxgridsize)=0.5d0*fr(ir)*auxgrid(ir)
!           do ir=2,auxgridsize
!             gr(auxgridsize)=0.5d0*(fr(ir)+fr(ir-1))*(auxgrid(ir)-auxgrid(ir-1))+gr(auxgridsize)
!           enddo
           Call fderiv (-1, auxgridsize, auxgrid, fr, gr, cf)
           ffacg (ig) = (fourpi/omega) * gr (auxgridsize)
         End Do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif

         deallocate(auxgrid,auxrho)

         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ig = 1, ngvec
               ifg = igfft (ig)
               zfft (ifg) = zfft (ifg) + ffacg (ig) * conjg (sfacg(ig, ias))
            End Do
         End Do
      End Do
      call timesec(tb)
      write(*,*) 'rhoinit, step 1:',tb-ta
      call timesec(ta)
! compute the tails in each muffin-tin

      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            zfmt (:, :) = 0.d0
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE) SHARED(zfft,is,ias,rcmt,gc,nrcmt,ngvec,sfacg,ffacg,ylmg,zfmt,igfft,vgc) PRIVATE(ig,irc,x,zt1,zt2,zt3,jj,z2fmt,boundlo,boundhi,ifg,yy,cs,sn,update,xi,jthr,r,tp,whichthread,nthreads)
#endif
            z2fmt (:, :) = 0.d0
#ifdef USEOMP
!$OMP DO
#endif
            Do ig = 1, ngvec
               ifg = igfft (ig)
!               Call sphcrd (vgc(:, ig), r, tp)
!               Call genylm (lmax, tp, yy)
!               yy=conjg(yy)
              yy=conjg(ylmg(1:4,ig))
               zt1 = fourpi * zfft (ifg) * sfacg (ig, ias)

!               call timesec(ta)
! Here we are supposed to compute the spherical Bessel functions for every r.
! The thing that we actually want to do can be described by 4 lines.
!               Do irc = 1, nrcmt (is)
!                 x = gc (ig) * rcmt (irc, is)
!                 Call sbessel (lmax, x, jj(:,irc))
!               enddo
! Unfortunately, this is slow.

! To make the code faster, we explicitly use assumptions:
! * lmax=1,
! * rcmt(:,is) is sorted.

! First, we separate x<1d-8 from large x.
! It is done using binary search.
               jthr=1d-8/gc (ig)
               boundhi=nrcmt (is)
               boundlo=1
               if (rcmt (1, is).gt.jthr) then
                 boundlo=0
                 boundhi=1
               elseif (rcmt (nrcmt (is), is).lt.jthr) then
                 boundlo=nrcmt (is)
                 boundhi=nrcmt (is)+1
               else
                 do while (boundhi-boundlo.gt.1)
                   irc=(boundhi+boundlo)/2
                   if (rcmt (irc, is).gt.jthr) then
                     boundhi=irc
                   else
                     boundlo=irc
                   endif
                 enddo
               endif
! Second, we apply the Taylor expansion for small x.
               do irc=1,boundlo
                 x=gc (ig)*rcmt (irc, is)
                 jj(0,irc)=1d0-sixth*x**2
                 jj(1,irc)=third*x*(1d0-0.1d0*x**2)
               enddo
! Third, we apply the actual formula for spherical Bessel functions for large x.
! j0(x)=sin(x)/x
! j1(x)=sin(x)/x**2-cos(x)/x
               do irc=boundhi,nrcmt(is)
                 x=gc (ig)*rcmt (irc, is)
                 xi=1d0/x
                 cs=cos(x)
                 sn=sin(x)
                 jj(0,irc)=sn*xi
                 jj(1,irc)=(jj(0,irc)-cs)*xi
               enddo
! End of story. :)

               update(1)=zt1*yy(1)
               update(2:4)=zt1*yy(2:4)*zi
               Do irc = 1, nrcmt (is)
                z2fmt (1, irc) = z2fmt (1, irc) + jj(0,irc) * update(1)
                z2fmt (2:4, irc) = z2fmt (2:4, irc) + jj(1,irc) * update(2:4)
               End Do
            End Do
#ifdef USEOMP
!$OMP END DO
            nthreads=omp_get_num_threads()
            whichthread=omp_get_thread_num()
            do irc=0,nthreads-1
              if (irc.eq.whichthread) then
                zfmt(1:4,1:nrcmt (is))= zfmt(1:4,1:nrcmt (is))+ z2fmt(1:4,1:nrcmt (is))
              endif
!$OMP BARRIER
            enddo
!$OMP END PARALLEL
#else
            zfmt=z2fmt
#endif
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
! Actually this is sloppy, since we have to subtract the contribution that is described by plane waves.
! Instead we subtract the model function that we try to expand in plane waves. In the limit of large Gmax,
! both things are equivalent. But should we even care? This is only the initial guess.
               t2 = (t1+sprho(ir, is)-a(is)*spr(ir,is)**2-c(is)) / y00 
               rhomt (1, ir, ias) = rhomt (1, ir, ias) + t2
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
      Deallocate (jlgr, ffacg, zfmt, zfft,z2fmt,a,c,jj)
      call timesec(tb)
      write(*,*) 'rhoinit, step 3:',tb-ta
!      stop
      Return
End Subroutine
!EOC
