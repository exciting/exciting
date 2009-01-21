
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine xsestimate
  use modmain
  use modxs
  use m_genfilname
  use m_getunit
  use m_gndstateq
  implicit none
  real(8) :: d_ev(2),d_apwcmt(2),d_pmat(2),d_emat(2),d_x0(2),d_tetw(2)
  real(8) :: d_fxcbse,d_sci,d_eci,d_tot(2)
  real(8) :: m_bseham,m_bseeigvec
  real(8) :: mb,gb
  integer :: sreal,scmplx,nkkpt,bsehamsiz,un
  call init0
  call init1
  call init2
  nst2=nempty+1
  nst1=nstsv-nst2
  ! estimate disk space usage
  sreal = 8
  scmplx = 16
  mb=1024**2
  gb=1024*mb
  nkkpt=(nkpt*(nkpt+1))/2
  bsehamsiz=nkpt*nbfbse*nafbse
  ! eigenvectors (spin-unpolarized)
  d_ev(1) = dble(scmplx*nkpt)*dble(nmatmax*nstfv)
  d_ev(2) = d_ev(1)*dble(nqpt+1)
  ! tetrahedron weights
  d_tetw(1) = dble(sreal*nkpt*nwdf)*dble(nst1*nst2*3)
  d_tetw(2) = d_tetw(1)
  ! expansion coefficients of muffin-tin wavefunctions
  d_apwcmt(1) = dble(scmplx*nkpt)*dble(nstfv*apwordmax*lmmaxapw*natmtot)
  d_apwcmt(2) = d_apwcmt(1)*dble(nqpt+1)
  ! matrix elements of the momentum operator
  d_pmat(1) = dble(scmplx*nkpt)*dble(3*nstsv*nstsv)
  d_pmat(2) = d_pmat(1)
  ! matrix elements of the plane wave
  d_emat(1) = dble(scmplx*nkpt)*dble(2*nst1*nst2*ngqmax)
  d_emat(2) = d_emat(1)*dble(nqpt)
  ! Kohn-Sham response function
  d_x0(1) = dble(scmplx*nwdf)*dble(ngqmax**2+3*2*ngqmax+3)
  d_x0(2) = d_x0(1)*dble(nqpt)
  ! BSE fxc-kernel
  d_fxcbse=d_x0(1)
  ! screened and exchange Coulomb interaction
  d_sci=scmplx*nkkpt*(nst1*nst2)**2
  d_eci=d_sci
  ! totals
  d_tot(1) = d_ev(1)+d_tetw(1)+d_apwcmt(1)+d_pmat(1)+d_emat(1)+d_x0(1)+d_sci+ &
       d_eci+d_fxcbse
  d_tot(2) = d_ev(2)+d_tetw(2)+d_apwcmt(2)+d_pmat(2)+d_emat(2)+d_x0(2)+d_sci+ &
       d_eci+d_fxcbse
  ! memory consumption
  m_bseham=scmplx*bsehamsiz**2
  m_bseeigvec=scmplx*bsehamsiz*nexcitmax
  !write information to file
  call getunit(un)
  open(un,file='SIZES.OUT',form='formatted',action='write',&
       status='replace')
  write(un,*)
  write(un,'(a)') 'Relevant parameters:'
  write(un,'(a,i6)') ' nwdf      :',nwdf
  write(un,'(a,i6)') ' nqpt      :',nqpt
  write(un,'(a,i6)') ' ngqmax    :',ngqmax
  write(un,'(a,i6)') ' nkpt      :',nkpt
  write(un,'(a,i6)') ' nkkpt     :',nkkpt
  write(un,'(a,i6)') ' nstsv     :',nstsv
  write(un,'(a,i6)') ' nst1      :',nst1
  write(un,'(a,i6)') ' nst2      :',nst2
  write(un,'(a,i6)') ' nbfbse    :',nbfbse
  write(un,'(a,i6)') ' nafbse    :',nafbse
  write(un,'(a,i6)') ' nmatmax   :',nmatmax
  write(un,'(a,i6)') ' ngkmax    :',ngkmax
  write(un,'(a,i6)') ' nlotot    :',nlotot
  write(un,'(a,i6)') ' apwordmax :',apwordmax
  write(un,'(a,i6)') ' lmmaxapw  :',lmmaxapw
  write(un,'(a,i6)') ' natmtot   :',natmtot
  write(un,*)
  write(un,'(a)') 'Other parameters:'  
  write(un,'(a,i8)') ' BSE Hamiltonian size  :',bsehamsiz
  write(un,*)
  write(un,'(a)') 'Estimated memory consumption in GB (one q-point/&
       &total):'
  write(un,'(a,2f12.3)') ' BSE Hamiltonian         :',m_bseham/gb
  write(un,'(a,2f12.3)') ' BSE eigenvectors        :',m_bseeigvec/gb
  write(un,'(a,2f12.3)') ' BSE totals *            :',(m_bseham+m_bseeigvec)/gb
  write(un,*)
  write(un,'(a)') 'Estimated disk usage for large files in GB (one q-point/&
       &total):'
  write(un,'(a,2f12.3)') ' eigenvectors                 :',d_ev(1)/gb,d_ev(2)/gb
  write(un,'(a,2f12.3)') ' tetrahedron weights          :',d_tetw(1)/gb, &
       d_tetw(2)/gb
  write(un,'(a,2f12.3)') ' MT expansion coeffs          :',d_apwcmt(1)/gb,&
       d_apwcmt(2)/gb
  write(un,'(a,2f12.3)') ' matr. el. of mom. op.        :',d_pmat(1)/gb, &
       d_pmat(2)/gb
  write(un,'(a,2f12.3)') ' matr. el. of plane wave      :',d_emat(1)/gb, &
       d_emat(2)/gb
  write(un,'(a,2f12.3)') ' KS response function         :',d_x0(1)/gb,d_x0(2)/gb
  write(un,'(a,2f12.3)') ' BSE-kernel                   :',d_fxcbse/gb
  write(un,'(a,2f12.3)') ' screened Coulomb interaction :',d_sci/gb
  write(un,'(a,2f12.3)') ' exchange Coulomb interaction :',d_eci/gb
  write(un,'(a,2f12.3)') '-----------------------------------------------------&
       &------------'
  write(un,'(a,2f12.3)') ' total                        :',d_tot(1)/gb, &
       d_tot(2)/gb
  write(un,*)
  close(un)
end subroutine xsestimate
