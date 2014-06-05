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
      Real (8) :: x, t1, t2, jlgr01, jlgr(0:1)
      Complex (8) zt1, zt2, zt3,yy(4),update(4)
      Real (8) :: ta,tb, tc,td
! automatic arrays
      Real (8) :: fr (spnrmax), gr (spnrmax), cf (3, spnrmax)
! allocatable arrays
      Real (8), Allocatable :: jj(:,:)
      Real (8), Allocatable :: th (:, :)
      Real (8), Allocatable :: ffacg (:)
      Real (8), Allocatable :: rhomodel(:,:)
      Complex (8), Allocatable :: zfmt (:, :),z2fmt (:,:)
      Complex (8), Allocatable :: zfft (:)

      integer :: auxgridsize,mtgridsize,lastpoint
      Real (8), Allocatable :: auxgrid(:),auxrho(:),a(:),c(:),b(:)
      Real (8) :: rhoder,rhoder2, tp(2),r
      real(8), parameter :: threshold=1d-12
      integer, parameter :: PointsPerPeriod=20

      integer:: boundhi,boundlo,middle
      real(8):: jthr,cs,sn,xi
#ifdef USEOMP
      integer:: whichthread,nthreads
#endif
      real(8), parameter :: sixth=1d0/6d0
      real(8), parameter :: third=1d0/3d0
!      real(8) :: cpu_seconds
!      external :: cpu_seconds
      integer :: nsw,isw,jsw,info  ! number of spherical waves
      real(8),allocatable :: swc(:),swoverlap(:,:),sine(:),cosine(:),swgr(:),pwswc(:,:),swoverlap2(:,:,:),pwswc2(:)
      complex(8),allocatable :: swc2(:,:,:),swctmp(:,:,:)!,pwswc(:,:)
      real(8) :: maxswg,swg,rhotest
      real(8) :: aa,bb,ans,rmt3


! maximum angular momentum for density initialisation
      lmax = 1
      lmmax = (lmax+1) ** 2
! allocate local arrays
!      Allocate (jj(0:lmax, nrcmtmax))
!      Allocate (jlgr(0:lmax))
      Allocate (ffacg(ngvec))
      Allocate (zfmt(lmmax, nrcmtmax))
      Allocate (rhomodel(nrcmtmax,nspecies))
!      Allocate (z2fmt(lmmax, nrcmtmax))
      Allocate (zfft(ngrtot))
      allocate(a(nspecies))
      allocate(b(nspecies))
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
! Do we still need this new grid? Maybe not so much...
         mtgridsize=PointsPerPeriod*int(rmt(is)*input%groundstate%gmaxvr/(2d0*pi))+1
         lastpoint=spnr(is)
         do while ((lastpoint.gt.nrmt(is)).and.(sprho (lastpoint, is).lt.threshold))
           lastpoint=lastpoint-1
         enddo

         auxgridsize=mtgridsize+lastpoint-nrmt(is)!+1
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
! pick .false. if you need extra smoothness 
if (.false.) then 
!-------- quadratic function
         a(is)=rhoder/(2d0*rmt(is))
         c(is)=sprho (nrmt(is), is)-a(is)*rmt(is)**2
! Interior of MT is filled with a*r**2+c
         do ir = 1, mtgridsize
           auxrho(ir)=a(is)*auxgrid(ir)**2+c(is)
         enddo
else
!-------- biquadratic function
         a(is)=(rhoder2*rmt(is)-rhoder)/(8d0*rmt(is)**3)
         b(is)=(3d0*rhoder-rhoder2*rmt(is))/(4d0*rmt(is))
         c(is)=sprho (nrmt(is), is)-a(is)*rmt(is)**4-b(is)*rmt(is)**2
! Interior of MT is filled with a*r**2+c
         do ir = 1, mtgridsize
           auxrho(ir)=a(is)*auxgrid(ir)**4+b(is)*auxgrid(ir)**2+c(is)
         enddo
endif

! Exterior is filled with the actual atomic density.
         auxrho(mtgridsize+1:auxgridsize)=sprho(nrmt(is)+1:lastpoint, is)

! Initialise auxiliary basis - spherical Bessel functions
         nsw=int(2d0*input%groundstate%gmaxvr*auxgrid(auxgridsize)/(pi))+1
         allocate(swc(0:nsw))
         allocate(swgr(0:nsw))
         allocate(swoverlap(0:nsw,0:nsw))
         maxswg=2d0*input%groundstate%gmaxvr ! 2*pi*dble(nsw)/auxgrid(auxgridsize) !input%groundstate%gmaxvr
         allocate(sine(0:nsw))
         allocate(cosine(0:nsw))


         do isw=0,nsw
           aa=maxswg*dble(isw)/dble(nsw)*auxgrid(auxgridsize)
           swgr(isw)=aa
           sine(isw)=sin(aa)
           cosine(isw)=cos(aa)
         enddo

! Compute overlaps of the auxiliary basis functions - spherical Bessel functions
         do isw=0,nsw
           aa=maxswg*dble(isw)/dble(nsw)*auxgrid(auxgridsize)
           do jsw=0,nsw
             bb=maxswg*dble(jsw)/dble(nsw)*auxgrid(auxgridsize)
             if (isw.lt.jsw) then
               if (isw.eq.0) then
                 swoverlap(isw,jsw)=auxgrid(auxgridsize)**3*(sine(jsw)-bb*cosine(jsw))/(bb**3)
                 swoverlap(jsw,isw)=swoverlap(isw,jsw)
               else
                 swoverlap(isw,jsw)=auxgrid(auxgridsize)**3*(bb*cosine(jsw)*sine(isw)-aa*cosine(isw)*sine(jsw))/(aa**3*bb-aa*bb**3)
                 swoverlap(jsw,isw)=swoverlap(isw,jsw)
               endif
             elseif (isw.eq.jsw) then
               if (isw.eq.0) then
                 swoverlap(isw,isw)=auxgrid(auxgridsize)**3/3d0
               else
                 swoverlap(isw,jsw)=auxgrid(auxgridsize)**3*(aa-cosine(isw)*sine(isw))/(2*aa**3)
                 swoverlap(jsw,isw)=swoverlap(isw,jsw)
               endif
             endif
           enddo
         enddo
! Factorisation is necessary for the least square fit in the L^2 sense
                  call dpotrf('U',&                       ! upper or lower part
                       nsw+1, &               ! size of matrix
                       swoverlap, &                 ! matrix
                       nsw+1, &               ! leading dimension
                       info &                     ! error message
                      )

! Density overlap with the auxilliary basis functions
         do ig=0,nsw
           swg=maxswg*dble(ig)/dble(nsw)
           do ir = 1,auxgridsize
             x=swg*auxgrid(ir)
             Call sbessel (0, x, jlgr01)
             fr (ir) = auxrho(ir) * jlgr01 * auxgrid(ir) ** 2
           enddo
           Call fderiv (-1, auxgridsize, auxgrid, fr, gr, cf)
           swc(ig)=gr(auxgridsize)
         enddo
! Density is expanded in the auxilliary basis.
! To do it, we find coefficients that minimise 
! \sum_i C_i \int_0^R j_0(k_i r) \rho(r) r^2 dr. 
            call dpotrs('U', &                      ! upper or lower part
                         nsw+1,  &                      ! size
                         1,  &                      ! number of right-hand sides
                         swoverlap, &      ! factorized matrix
                         nsw+1, &                       ! leading dimension
                         swc(0), &                    ! right-hand side / solution
                         nsw+1, &                       ! leading dimension
                         info &                     ! error message
                       )
! Calculate the smooth model density. 
! If the number of basis fns is large enough it should coincide with a*r^2+c.

        do irc=1,nrcmt(is)
          rhotest=0d0
          do ig=0,nsw
            x=maxswg*dble(ig)/dble(nsw)*rcmt (irc, is)
            call sbessel (0, x, jlgr01)
            rhotest=rhotest+swc(ig)*jlgr01
          enddo
          rhomodel(irc,is)=rhotest
        enddo

! Expand the model density in plane waves.
! To reduce the effort, we express the model density in terms of the auxilliary basis.
        ffacg=0d0
        ffacg(1)=ffacg(1)+third*swc(0)
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(none) PRIVATE(ig,bb) SHARED(ngvec,auxgrid,auxgridsize,gc,ffacg,swc)
!$OMP DO 
#endif
        do ig=2,ngvec
          bb=gc(ig)*auxgrid(auxgridsize)
          ffacg(ig)=ffacg(ig)+(sin(bb)-bb*cos(bb))/(bb**3)*swc(0)
        enddo
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif

        do isw=1,nsw
          swg=maxswg*dble(isw)/dble(nsw)
          aa=swg*auxgrid(auxgridsize)
          ffacg(1)=ffacg(1)+(sine(isw)-aa*cosine(isw))/(aa**3)*swc(isw)
        enddo

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(none) PRIVATE(ig,bb,sn,cs,aa,isw) SHARED(ngvec,auxgrid,auxgridsize,swgr,gc,ffacg,sine,cosine,swc,nsw)
!$OMP DO 
#endif
          do ig=2,ngvec
            bb=gc(ig)*auxgrid(auxgridsize)
            sn=sin(bb)
            cs=cos(bb)
            if (abs(dnint (bb/swgr(1))-bb/swgr(1)).lt.1d-8) then
              do isw=1,nsw
               aa=swgr(isw)
               if (abs(bb-aa).lt.1d-8) then
                 ffacg(ig)=ffacg(ig)+(aa-cosine(isw)*sine(isw))/(2*aa**3)*swc(isw)
               else
                 ffacg(ig)=ffacg(ig)+(bb*cs*sine(isw)-aa*cosine(isw)*sn)/(aa**3*bb-aa*bb**3)*swc(isw)
               endif
             enddo
            else
             do isw=1,nsw
               aa=swgr(isw)
               ffacg(ig)=ffacg(ig)+(bb*cs*sine(isw)-aa*cosine(isw)*sn)/(aa**3*bb-aa*bb**3)*swc(isw)
             enddo
            endif
          enddo
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif

        deallocate(cosine,sine,swgr)
        deallocate(swc,swoverlap)
        ffacg=ffacg*(fourpi/omega)*auxgrid(auxgridsize)**3
        deallocate(auxgrid,auxrho)

! Sum all the contribution from different atoms that correspond to same species.
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ig = 1, ngvec
               ifg = igfft (ig)
               zfft (ifg) = zfft (ifg) + ffacg (ig) * conjg (sfacg(ig, ias))
            End Do
         End Do
      End Do
      call timesec(tb)
      !write(*,*) 'rhoinit, step 1:',tb-ta
      call timesec(ta)
! compute the tails in each muffin-tin

! Choose .false. if you want to revert to old ways
if (.true.) then
      Do is = 1, nspecies
         nsw=int(2*input%groundstate%gmaxvr*rmt(is)/(pi))+1
         allocate(swc2(0:nsw,4,natoms(is)))
!         allocate(pwswc(0:nsw,2))
!         allocate(pwswc2(0:nsw))
         allocate(swgr(0:nsw))
         allocate(swoverlap2(0:nsw,0:nsw,2))
         maxswg=2d0*input%groundstate%gmaxvr ! 2*pi*dble(nsw)/auxgrid(auxgridsize) !input%groundstate%gmaxvr
         allocate(sine(0:nsw))
         allocate(cosine(0:nsw))


         do isw=0,nsw
           aa=maxswg*dble(isw)/dble(nsw)*rmt(is)
           swgr(isw)=aa
           sine(isw)=sin(aa)
           cosine(isw)=cos(aa)
         enddo
         rmt3=rmt(is)**3

! Compute overlaps of the auxiliary basis functions - spherical Bessel functions up to l=1
         do isw=0,nsw
           aa=maxswg*dble(isw)/dble(nsw)*rmt(is)
           do jsw=0,nsw
             bb=maxswg*dble(jsw)/dble(nsw)*rmt(is)
             if (isw.lt.jsw) then
               if (isw.eq.0) then
                 swoverlap2(isw,jsw,1)=rmt3*(sine(jsw)-bb*cosine(jsw))/(bb**3)
                 swoverlap2(jsw,isw,1)=swoverlap2(isw,jsw,1)

                 swoverlap2(isw,jsw,2)=0d0
                 swoverlap2(jsw,isw,2)=0d0

               else
                 swoverlap2(isw,jsw,1)=rmt3*(bb*cosine(jsw)*sine(isw)-aa*cosine(isw)*sine(jsw))/(aa**3*bb-aa*bb**3)
                 swoverlap2(jsw,isw,1)=swoverlap2(isw,jsw,1)

                 swoverlap2(isw,jsw,2)=rmt3*(aa**2*bb*cosine(jsw)*sine(isw)    &
                                                  -aa*bb**2*cosine(isw)*sine(jsw)    &
                                                  -(aa**2-bb**2)*sine(isw)*sine(jsw)) &
                                                  /(aa**2*bb**2*(aa**2-bb**2))
                 swoverlap2(jsw,isw,2)=swoverlap2(isw,jsw,2)
               endif
             elseif (isw.eq.jsw) then
               if (isw.eq.0) then
                 swoverlap2(isw,isw,1)=rmt3/3d0
                 swoverlap2(isw,jsw,2)=0d0
               else
                 swoverlap2(isw,jsw,1)=rmt3*(aa-cosine(isw)*sine(isw))/(2*aa**3)
                 swoverlap2(jsw,isw,1)=swoverlap2(isw,jsw,1)
                 swoverlap2(isw,jsw,2)=rmt3*(aa**2-2d0*sine(isw)**2 +aa*cosine(isw)*sine(isw))/(2*aa**4)
                 swoverlap2(jsw,isw,2)=swoverlap2(isw,jsw,2)
               endif
             endif
           enddo
         enddo


! Invert overlap matrix - prepare for the least square fitting
         swoverlap2(0,:,2)=0d0
         swoverlap2(:,0,2)=0d0
         swoverlap2(0,0,2)=1d0
         call dpotrf('U',&                       ! upper or lower part
                     nsw+1, &               ! size of matrix
                     swoverlap2(0,0,1), &                 ! matrix
                     nsw+1, &               ! leading dimension
                     info &                     ! error message
                    )
         call dpotri('U',&                       ! upper or lower part
                     nsw+1, &               ! size of matrix
                     swoverlap2(0,0,1), &                 ! matrix
                     nsw+1, &               ! leading dimension
                     info &                     ! error message
                    )
         do isw=0,nsw
           swoverlap2(:,isw,1)=swoverlap2(isw,:,1)
         enddo
         call dpotrf('U',&                       ! upper or lower part
                     nsw+1, &               ! size of matrix
                     swoverlap2(0,0,2), &                 ! matrix
                     nsw+1, &               ! leading dimension
                     info &                     ! error message
                    )
         call dpotri('U',&                       ! upper or lower part
                     nsw+1, &               ! size of matrix
                     swoverlap2(0,0,2), &                 ! matrix
                     nsw+1, &               ! leading dimension
                     info &                     ! error message
                    )
         swoverlap2(0,0,2)=0d0
         do isw=0,nsw
           swoverlap2(:,isw,2)=swoverlap2(isw,:,2)
         enddo


! Express plane waves in terms of spherical Bessel functions
! and calculate them of the radial grid

        swc2(:,:,:)=zzero
        ifg = igfft (1)

! Consider G=0 separately
        yy=conjg(ylmg(1:4,1))
        do ia= 1, natoms (is)
         ias = idxas (ia, is)
         swc2(0,1,ia)=zfft (ifg)*yy(1)*fourpi* sfacg (1, ias)
        enddo

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(none) PRIVATE(ig,bb,sn,cs,aa,isw,yy,ifg,pwswc,pwswc2,swctmp,ia,ias,irc,whichthread,nthreads) SHARED(ngvec,swgr,gc,sine,cosine,nsw,ylmg,natoms,sfacg,zfft,swc2,rmt,is,igfft,swoverlap2,idxas,rmt3)
        allocate(swctmp(0:nsw,4,natoms(is)))
        swctmp=zzero
#endif
        allocate(pwswc(0:nsw,2))
        allocate(pwswc2(0:nsw))
#ifdef USEOMP
!$OMP DO 
#endif
! Now other Gs follow
        do ig=2,ngvec
          ifg = igfft (ig)
          yy=conjg(ylmg(1:4,ig))
          bb=gc(ig)*rmt(is)
          sn=sin(bb)
          cs=cos(bb)
! Plane wave expands into spherical Bessel functions, 
! but we don't want to compute them for every G on each 
! radial grid point.
! Insted we expand each occuring spherical Bessel function 
! in terms of the auxilliary basis (also spherical Bessel 
! functions) with the size that is much smaller than the number Gs.

! First, we compute overlaps
          pwswc(0,1)=rmt3*(sn-bb*cs)/(bb**3)
          pwswc(0,2)=0d0
          if (abs(dnint(bb/swgr(1))-bb/swgr(1)).lt.1d-8) then
            do isw=1,nsw
              aa=swgr(isw)
              if (abs(bb-aa).lt.1d-8) then
                pwswc(isw,1)=rmt3*(aa-cosine(isw)*sine(isw))/(2*aa**3)
                pwswc(isw,2)=rmt3*(aa**2-2d0*sine(isw)**2 +aa*cosine(isw)*sine(isw))/(2*aa**4)
              else
                pwswc(isw,1)=rmt3*(bb*cs*sine(isw)-aa*cosine(isw)*sn)/(aa**3*bb-aa*bb**3)
                pwswc(isw,2)=rmt3*(aa**2*bb*cs*sine(isw)    &
                                   -aa*bb**2*cosine(isw)*sn    &
                                   -(aa**2-bb**2)*sine(isw)*sn) &
                                   /(aa**2*bb**2*(aa**2-bb**2))
              endif
            enddo
          else
            do isw=1,nsw
              aa=swgr(isw)
              pwswc(isw,1)=rmt3*(bb*cs*sine(isw)-aa*cosine(isw)*sn)/(aa**3*bb-aa*bb**3)
              pwswc(isw,2)=rmt3*(aa**2*bb*cs*sine(isw)    &
                                -aa*bb**2*cosine(isw)*sn    &
                                -(aa**2-bb**2)*sine(isw)*sn) &
                                /(aa**2*bb**2*(aa**2-bb**2))
            enddo
          endif

! Second, we use the overlaps to obtain the least-square coefficients
          call dgemv('N',nsw+1,nsw+1,1d0,swoverlap2(0,0,1),nsw+1,pwswc(0,1),1,0d0,pwswc2,1)
          pwswc(:,1)=pwswc2
          call dgemv('N',nsw+1,nsw+1,1d0,swoverlap2(0,0,2),nsw+1,pwswc(0,2),1,0d0,pwswc2,1)
          pwswc(:,2)=pwswc2
! pwswc contains the expansion of j_n(G*r) in terms of the auxilliary basis

! Update the contribution of rho_G for every atom of the given species.
! This contribution is still expressed in terms of the auxilliary basis.
          do ia= 1, natoms (is)
            ias = idxas (ia, is)
#ifdef USEOMP
            swctmp(:,1,ia)=swctmp(:,1,ia)+pwswc(:,1)*zfft (ifg)*yy(1)*   fourpi* sfacg (ig, ias)
            swctmp(:,2,ia)=swctmp(:,2,ia)+pwswc(:,2)*zfft (ifg)*yy(2)*zi*fourpi* sfacg (ig, ias)
            swctmp(:,3,ia)=swctmp(:,3,ia)+pwswc(:,2)*zfft (ifg)*yy(3)*zi*fourpi* sfacg (ig, ias)
            swctmp(:,4,ia)=swctmp(:,4,ia)+pwswc(:,2)*zfft (ifg)*yy(4)*zi*fourpi* sfacg (ig, ias)
#else
            swc2(:,1,ia)=swc2(:,1,ia)+pwswc(:,1)*zfft (ifg)*yy(1)*   fourpi* sfacg (ig, ias)
            swc2(:,2,ia)=swc2(:,2,ia)+pwswc(:,2)*zfft (ifg)*yy(2)*zi*fourpi* sfacg (ig, ias)
            swc2(:,3,ia)=swc2(:,3,ia)+pwswc(:,2)*zfft (ifg)*yy(3)*zi*fourpi* sfacg (ig, ias)
            swc2(:,4,ia)=swc2(:,4,ia)+pwswc(:,2)*zfft (ifg)*yy(4)*zi*fourpi* sfacg (ig, ias)
#endif
          enddo


       enddo
#ifdef USEOMP
!$OMP END DO
            nthreads=omp_get_num_threads()
            whichthread=omp_get_thread_num()
            do irc=0,nthreads-1
              if (irc.eq.whichthread) then
                swc2=swctmp+swc2
              endif
!$OMP BARRIER
            enddo
            deallocate(swctmp)
#endif
            deallocate(pwswc,pwswc2)
#ifdef USEOMP
!$OMP END PARALLEL
#endif

! Calculate the density at each radial grid point on a coarse mesh
! and transform it to the fine mesh.
       do ia= 1, natoms (is)
         ias = idxas (ia, is)
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(none) PRIVATE(irc,yy,x,jlgr,isw) SHARED(nsw,nrcmt,ia,swc2,maxswg,rcmt,zfmt,is)
!$OMP DO 
#endif

         do irc=1,nrcmt(is)
           yy=0d0
           yy(1)=swc2(0,1,ia)
           do isw=1,nsw
             x=maxswg*dble(isw)/dble(nsw)*rcmt (irc, is)
             call sbessel (1, x, jlgr)
             yy(1)=yy(1)+swc2(isw,1,ia)*jlgr(0)
             yy(2:4)=yy(2:4)+swc2(isw,2:4,ia)*jlgr(1)
           enddo
           zfmt(1:4,irc)=yy(1:4)
         enddo
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif

         irc = 0
         Do ir = 1, nrmt (is), input%groundstate%lradstep
           irc = irc + 1
           Call ztorflm (lmax, zfmt(:, irc), rhomt(:, ir, ias))
         End Do
       enddo

       deallocate(cosine,sine,swgr,swc2,swoverlap2)
         

     enddo

else

      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            zfmt (:, :) = 0.d0
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE) SHARED(zfft,is,ias,rcmt,gc,nrcmt,ngvec,sfacg,ffacg,ylmg,zfmt,igfft,vgc) PRIVATE(ig,irc,x,zt1,zt2,zt3,jj,z2fmt,boundlo,boundhi,ifg,yy,cs,sn,update,xi,jthr,r,tp,whichthread,nthreads)
#endif
            Allocate (jj(0:1, nrcmt (is)))
            Allocate (z2fmt(4, nrcmt (is)))
            z2fmt (:, :) = 0.d0
#ifdef USEOMP
!$OMP DO
#endif
            Do ig = 1, ngvec
               ifg = igfft (ig)
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
            deallocate(jj,z2fmt)
!$OMP END PARALLEL
#else
            zfmt=z2fmt
            deallocate(jj,z2fmt)
#endif
            irc = 0
            Do ir = 1, nrmt (is), input%groundstate%lradstep
               irc = irc + 1
               Call ztorflm (lmax, zfmt(:, irc), rhomt(:, ir, ias))
            End Do
         End Do
      End Do

endif




      call timesec(tb)
      !write(*,*) 'rhoinit, step 2:',tb-ta
      call timesec(ta)
 
! Remove model density from rhomt
      Do is = 1, nspecies
        Do ia = 1, natoms (is)
          ias = idxas (ia, is)
          rhomt(1,1:nrcmt(is),ias)=rhomt(1,1:nrcmt(is),ias)-rhomodel(1:nrcmt(is),is)/y00
        enddo
      enddo

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
!               t2 = (t1+sprho(ir, is)-a(is)*spr(ir,is)**2-c(is)) / y00 
               t2 = (t1+sprho(ir, is)) / y00
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
      Deallocate (ffacg, zfmt, zfft,a,b,c,rhomodel)
      call timesec(tb)
      !write(*,*) 'rhoinit, step 3:',tb-ta
!      stop
      Return
End Subroutine
!EOC
