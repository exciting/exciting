! Copyright (C) 2007-2010 D. Nabok, P. Puschnig and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

module param

!  parameters
   real*8, parameter   :: eps  = 1.0d-10
   real*8, parameter   :: pi   = 3.1415926535d0
   real*8, parameter   :: bohr = 0.529177d0

!  density units
   character*180       :: phifile, xsffile, xsfgradients
   character*4         :: densunits

!  integration parameters
   integer    seed
   integer    flags
   integer*8  mineval, maxeval
   double precision epsrel, epsabs
   double precision userdata
   integer    key1, key2, key3, maxpass
   double precision border, maxchisq, mindeviation
   integer    ldxgiven 
   integer*8  ngiven, nextra
   
   parameter (seed = 0)
   parameter (userdata = 0.0d0)
!  parameter (epsrel = 1.D-01)
!  parameter (epsabs = 1.D-03)
   parameter (flags = 0)
   parameter (mineval = 96)
!  parameter (maxeval = 9223372036854775807)
!  parameter (key1 = 1000)
   parameter (key2 = 1)
   parameter (key3 = 1)
   parameter (maxpass = 5)
   parameter (border = 0.D0)
   parameter (maxchisq = 1.D0)
   parameter (mindeviation = .25D0)
   parameter (ngiven = 0)
   parameter (ldxgiven = 0)
   parameter (nextra = 0)
   
!  structure   
   real*8              :: origin(3)
   real*8              :: vectors(3,3)
   real*8              :: volume
   integer             :: nrx,nry,nrz

!  electron density and related quantities on a real space grid    
   integer             :: nx,ny,nz
   integer             :: ixstep,iystep,izstep   
   real*8, allocatable :: density(:,:,:)
   real*8, allocatable :: graddensity(:,:,:)  ! gradient of the density

   character*10        :: vdWDF_version
   real*8              :: Zab_04, Zab_10, Cfac
   parameter (Zab_04=-0.8491d0)
   parameter (Zab_10=-1.8870d0 )
   parameter (Cfac=0.0089d0)
 
!  variables related to the vdW kernel phi  
   integer             :: nd, ndelta
   real*8              :: dmin,dmax,dstep,Dcutoff
   real*8              :: deltamin,deltamax,deltastep
   real*8, allocatable :: kernel(:,:)

   real*8              :: const

   contains 
   
   subroutine delete_arrays
     deallocate(density)
     deallocate(graddensity)
     deallocate(kernel)
   end subroutine delete_arrays   
         
end module param
