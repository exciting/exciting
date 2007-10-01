
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phonon
use modmain
implicit none
! local variables
integer iq,is,ia,ip,js,ja,jas,jp,ka
integer nph,iph,i1,i2,i3,iv(3)
integer natoms0(maxspecies)
real(8) v1(3),v2(3),v3(3),t1,t2
real(8) deltaph0,avec0(3,3)
real(8) atposc0(3,maxatoms,maxspecies)
real(8) ftp(3,maxatoms,maxspecies)
real(8) df(3,maxatoms,maxspecies,0:1)
complex(8) zsum,zt1
! external functions
real(8) r3taxi
external r3taxi
!------------------------!
!     initialisation     !
!------------------------!
! require forces
tforce=.true.
! no primitive cell determination
primcell=.false.
! initialise universal variables
call init0
! initialise q-point dependent variables
call init2
! switch off automatic determination of muffin-tin radii
autormt=.false.
! no shifting of atomic basis allowed
tshift=.false.
! determine k-point grid size from maximum de Broglie wavelength
autokpt=.true.
! store input values
natoms0(1:nspecies)=natoms(1:nspecies)
deltaph0=deltaph
avec0(:,:)=avec(:,:)
atposc0(:,:,:)=0.d0
do is=1,nspecies
  do ia=1,natoms(is)
    atposc0(:,ia,is)=atposc(:,ia,is)
  end do
end do
!---------------------------------------!
!     compute dynamical matrix rows     !
!---------------------------------------!
10 continue
natoms(1:nspecies)=natoms0(1:nspecies)
! find a dynamical matrix to calculate
call dyntask(80,iq,is,ia)
! check to see if mass is considered infinite
if (spmass(is).le.0.d0) then
  do ip=1,3
    write(80,*)
    do js=1,nspecies
      do ja=1,natoms0(js)
        do jp=1,3
          write(80,'(2G18.10," : is = ",I4,", ia = ",I4,", ip = ",I4)') 0.d0, &
           0.d0,js,ja,jp
        end do
      end do
    end do
  end do
  close(80)
  goto 10
end if
if (sqrt(vql(1,iq)**2+vql(2,iq)**2+vql(3,iq)**2).lt.epslat) then
  nph=0
else
  nph=1
end if
! start from atomic densities
task=0
! begin loop over polarisations
do ip=1,3
! loop over phases (cos and sin displacements)
  do iph=0,nph
! restore input values
    natoms(1:nspecies)=natoms0(1:nspecies)
    avec(:,:)=avec0(:,:)
    atposc(:,:,:)=atposc0(:,:,:)
    deltaph=deltaph0
! generate the supercell
    call phcell(iph,iq,is,ia,ip)
! run the ground-state calculation
    call gndstate
! store the total force for the first displacement
    do js=1,nspecies
      do ja=1,natoms(js)
        jas=idxas(ja,js)
        ftp(:,ja,js)=forcetot(:,jas)
      end do
    end do
! restore input values
    natoms(1:nspecies)=natoms0(1:nspecies)
    avec(:,:)=avec0(:,:)
    atposc(:,:,:)=atposc0(:,:,:)
! double the displacement
    deltaph=deltaph0+deltaph0
! generate the supercell again
    call phcell(iph,iq,is,ia,ip)
! run the ground-state calculation again starting from previous density
    task=1
    call gndstate
! store force derivative for current phase
    do js=1,nspecies
      do ja=1,natoms(js)
        jas=idxas(ja,js)
        df(:,ja,js,iph)=-(forcetot(:,jas)-ftp(:,ja,js))/deltaph0
      end do
    end do
  end do
!---------------------------------------------!
!     Fourier transform force derivatives     !
!---------------------------------------------!
  write(80,*)
! restore input values
  natoms(1:nspecies)=natoms0(1:nspecies)
  avec(:,:)=avec0(:,:)
  atposc(:,:,:)=atposc0(:,:,:)
! zero displacement
  deltaph=0.d0
  call phcell(0,iq,is,ia,ip)
  t1=1.d0/dble(ngridq(1)*ngridq(2)*ngridq(3))
! begin loops over species, atoms and polarisations
  do js=1,nspecies
    do ja=1,natoms0(js)
      do jp=1,3
! Fourier transform the force derivatives
        zsum=0.d0
! loop over phases (cos and sin displacements)
        do iph=0,nph
! loop over R-vectors
          do i3=ngridq(3)/2-ngridq(3)+1,ngridq(3)/2
            do i2=ngridq(2)/2-ngridq(2)+1,ngridq(2)/2
              do i1=ngridq(1)/2-ngridq(1)+1,ngridq(1)/2
                v1(:)=dble(i1)*avec(:,1)+dble(i2)*avec(:,2)+dble(i3)*avec(:,3)
                t2=-(vqc(1,iq)*v1(1)+vqc(2,iq)*v1(2)+vqc(3,iq)*v1(3))
                zt1=cmplx(cos(t2),sin(t2),8)
                if (iph.eq.1) zt1=zt1*zi
                v2(:)=v1(:)+atposc0(:,ja,js)
! convert atomic position to supercell lattice coordinates
                call r3mv(ainv,v2,v2)
                call r3frac(epslat,v2,iv)
! locate atom in current supercell
                do ka=1,natoms(js)
                  v3(:)=atposl(:,ka,js)
                  call r3frac(epslat,v3,iv)
                  if (r3taxi(v2,v3).lt.epslat) then
                    zsum=zsum+zt1*df(jp,ka,js,iph)
                    goto 30
                  end if
                end do
30 continue
              end do
            end do
          end do
        end do
        zsum=t1*zsum
! write dynamical matrix row to file
        write(80,'(2G18.10," : is = ",I4,", ia = ",I4,", ip = ",I4)') zsum,js, &
         ja,jp
      end do
    end do
  end do
  call flushifc(80)
! end loop over polarisations
end do
close(80)
! delete the eigenvector files
call delevec
goto 10
end subroutine

