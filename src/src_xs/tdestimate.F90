
subroutine tdestimate
  use modmain
  use modtddft
  use m_getunit
  implicit none
  real(8) :: d_ev(2), d_apwdlm(2), d_pmat(2), d_emat(2), d_x0(2), d_tetw(2)
  real(8) :: d_tot(2), mb,gb
  integer :: sreal,scmplx,un
  
  call init0
  call init1
  call init2xs

  ! estimate disk space usage
  sreal = 8
  scmplx = 16
  mb=1024**2
  gb=1024*mb
  ! eigenvectors (spin-unpolarized)
  d_ev(1) = dble(scmplx*nkpt)*dble(nmatmax*nstfv)
  d_ev(2) = d_ev(1)*dble(nqpt+1)
  ! tetrahedron weights
  d_tetw(1) = dble(sreal*nkpt*nwdf)*dble(nstval*nstcon*3)
  d_tetw(2) = d_tetw(1)
  ! expansion coefficients of muffin-tin wavefunctions
  d_apwdlm(1) = dble(scmplx*nkpt)*dble(nstfv*apwordmax*lmmaxapw*natmtot)
  d_apwdlm(2) = d_apwdlm(1)*dble(nqpt+1)
  ! matrix elements of the momentum operator
  d_pmat(1) = dble(scmplx*nkpt)*dble(3*nstsv*nstsv)
  d_pmat(2) = d_pmat(1)
  ! matrix elements of the plane wave
  d_emat(1) = dble(scmplx*nkpt)*dble(2*nstval*nstcon*ngqmax)
  d_emat(2) = d_emat(1)*dble(nqpt)
  ! Kohn-Sham response function
  d_x0(1) = dble(scmplx*nwdf)*dble(ngqmax**2+3*2*ngqmax+3)
  d_x0(2) = d_x0(1)*dble(nqpt)

  ! totals
  d_tot(1) = d_ev(1)+d_tetw(1)+d_apwdlm(1)+d_pmat(1)+d_emat(1)+d_x0(1)
  d_tot(2) = d_ev(2)+d_tetw(2)+d_apwdlm(2)+d_pmat(2)+d_emat(2)+d_x0(2)

  !write information to file
  call getunit(un)
  open(un,file='TDESTIMATE.OUT',form='formatted',action='write',&
       status='replace')
  write(un,*)
  write(un,'(a)') 'Relevant parameters:'
  write(un,'(a,i6)') ' nwdf     :',nwdf
  write(un,'(a,i6)') ' nqpt     :',nqpt
  write(un,'(a,i6)') ' ngqmax   :',ngqmax
  write(un,'(a,i6)') ' nkpt     :',nkpt
  write(un,'(a,i6)') ' nstsv    :',nstsv
  write(un,'(a,i6)') ' nstval   :',nstval
  write(un,'(a,i6)') ' nstcon   :',nstcon
  write(un,'(a,i6)') ' nmatmax  :',nmatmax
  write(un,'(a,i6)') ' apwordmax:',apwordmax
  write(un,'(a,i6)') ' lmmaxapw :',lmmaxapw
  write(un,'(a,i6)') ' natmtot  :',natmtot
  write(un,*)
  write(un,'(a)') 'Estimated disk usage for large files in GB (one q-point/&
       &total):'
  write(un,'(a,2f12.3)') ' eigenvectors           :',d_ev(1)/gb,d_ev(2)/gb
  write(un,'(a,2f12.3)') ' tetrahedron weights    :',d_tetw(1)/gb,d_tetw(2)/gb
  write(un,'(a,2f12.3)') ' MT expansion coeffs    :',d_apwdlm(1)/gb,&
       d_apwdlm(2)/gb
  write(un,'(a,2f12.3)') ' matr. el. of mom. op.  :',d_pmat(1)/gb,d_pmat(2)/gb
  write(un,'(a,2f12.3)') ' matr. el. of plane wave:',d_emat(1)/gb,d_emat(2)/gb
  write(un,'(a,2f12.3)') ' KS response function   :',d_x0(1)/gb,d_x0(2)/gb
  write(un,'(a,2f12.3)') '----------------------------------------------------------------'
  write(un,'(a,2f12.3)') ' total                  :',d_tot(1)/gb,d_tot(2)/gb
  write(un,*)
  close(un)

end subroutine tdestimate
