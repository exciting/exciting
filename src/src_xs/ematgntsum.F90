


! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine ematgntsum(iq, igq)
  use modmain
use modinput
  use modxs
  use m_findgntn0
  use m_getunit
  implicit none
  ! arguments
  integer, intent(in) :: iq, igq
  ! local variables
  integer::is, ia, ias
  integer::l1, l2, l3, m2, lm2
  integer::ilo, ilo1, ilo2, io, io1, io2
  integer::lmax1, lmax2, lmax3, lmmax1, lmmax2, lmmax3
  integer::u1, u2, u3, u4
  integer::m1, m3, lm1, lm3, cl1, cm1, cl2, cm2, cl3, cm3
  lmax1=max(input%xs%lmaxapwwf, lolmax)
  lmax2=input%xs%lmaxemat
  ! lmax1 and lmax3 should be the same
  lmax3=lmax1
  lmmax1=(lmax1+1)**2
  lmmax2=(lmax2+1)**2
  lmmax3=(lmax3+1)**2
  ! allocate arrays for radial integrals and Bessel functions
  if (allocated(intrgaa)) deallocate(intrgaa)
  if (allocated(intrgloa)) deallocate(intrgloa)
  if (allocated(intrgalo)) deallocate(intrgalo)
  if (allocated(intrglolo)) deallocate(intrglolo)
  allocate(intrgaa(lmmax1, apwordmax, lmmax3, apwordmax, natmtot))
  allocate(intrgloa(-lolmax:lolmax, nlomax, lmmax3, apwordmax, natmtot))
  allocate(intrgalo(-lolmax:lolmax, nlomax, lmmax3, apwordmax, natmtot))
  allocate(intrglolo(-lolmax:lolmax, nlomax, -lolmax:lolmax, nlomax, natmtot))
  ! allocate temporary arrays
  intrgaa(:, :, :, :, :)=zzero
  intrgloa(:, :, :, :, :)=zzero
  intrgalo(:, :, :, :, :)=zzero
  intrglolo(:, :, :, :, :)=zzero
  if (input%xs%dbglev.gt.2) then
     ! APW-APW
     call getunit(u1)
     open(unit = u1, file = 'IRADGAUNTaa'//filext, form = 'formatted', &
	  action = 'write', status = 'replace')
     write(u1, '(a)') 'igq, ias, lm1, io1, lm3, io2,   intrgaa'
     write(u1, '(a)') '------------------------------------------------------'
     ! APW-lo
     call getunit(u2)
     open(unit = u2, file = 'IRADGAUNTalo'//filext, form = 'formatted', &
	  action = 'write', status = 'replace')
     write(u2, '(a)') 'igq, ias, m3, ilo, lm1, io,     intrgalo'
     write(u2, '(a)') '------------------------------------------------------'
     ! lo-APW
     call getunit(u3)
     open(unit = u3, file = 'IRADGAUNTloa'//filext, form = 'formatted', &
	  action = 'write', status = 'replace')
     write(u3, '(a)') 'igq, ias, m1, ilo, lm3, io,     intrgloa'
     write(u3, '(a)') '------------------------------------------------------'
     ! lo-lo
     call getunit(u4)
     open(unit = u4, file = 'IRADGAUNTlolo'//filext, form = 'formatted', &
	  action = 'write', status = 'replace')
     write(u4, '(a)') 'igq, ias, m1, ilo1, m3, ilo2,   intrglolo'
     write(u4, '(a)') '------------------------------------------------------'
  end if
  ! begin loop over species
  do is=1, nspecies
     ! begin loop over atoms
     do ia=1, natoms(is)
	ias=idxas(ia, is)
        !---------------------------!
        !     APW-APW integrals     !
        !---------------------------!
	do cl1=1, l1shape
	   l1=l1map(cl1)
	   do cm1=1, m1shape(l1)
	      m1=m1map(l1, cm1)
	      lm1=idxlm(l1, m1)
	      do io1=1, apword(l1, is)
		 do cl2=1, l2shape(l1, m1)
		    l3=l2map(l1, m1, cl2)
		    do cm2=1, m2shape(l1, m1, l3)
		       m3=m2map(l1, m1, l3, cm2)
		       lm3=idxlm(l3, m3)
		       do io2=1, apword(l3, is)
			  do cl3=1, l3shape(l1, m1, l3, m3)
			     l2=l3map(l1, m1, l3, m3, cl3)
			     do cm3=1, m3shape(l1, m1, l3, m3, l2)
				m2=m3map(l1, m1, l3, m3, l2, cm3)
				lm2=idxlm(l2, m2)
				intrgaa(lm1, io1, lm3, io2, ias)= &
				     intrgaa(lm1, io1, lm3, io2, ias) &
				     +conjg(zil(l2))* &
				     riaa(l1, io1, l3, io2, l2, ias, igq)* &
				     conjg(ylmgq(lm2, igq, iq))* &
				     xsgnt(lm1, lm2, lm3)
			     end do
			  end do
			  if (input%xs%dbglev.gt.2) then
			     write(u1, '(6i5, 2g18.10)') &
				  igq, ias, lm1, io1, lm3, io2, &
				  intrgaa(lm1, io1, lm3, io2, ias)
			  end if
		       end do
		    end do
		 end do
	      end do
	   end do
	end do
        !-------------------------------------!
        !     APW-local-orbital integrals     !
        !-------------------------------------!
	do cl1=1, l1shape
	   l1=l1map(cl1)
	   do cm1=1, m1shape(l1)
	      m1=m1map(l1, cm1)
	      lm1=idxlm(l1, m1)
	      do io=1, apword(l1, is)
		 do ilo=1, nlorb(is)
		    l3=lorbl(ilo, is)
		    do cm2=1, m2shape(l1, m1, l3)
		       m3=m2map(l1, m1, l3, cm2)
		       lm3=idxlm(l3, m3)
		       do cl3=1, l3shape(l1, m1, l3, m3)
			  l2=l3map(l1, m1, l3, m3, cl3)
			  do cm3=1, m3shape(l1, m1, l3, m3, l2)
			     m2=m3map(l1, m1, l3, m3, l2, cm3)
			     lm2=idxlm(l2, m2)
			     intrgalo(m3, ilo, lm1, io, ias)= &
				  intrgalo(m3, ilo, lm1, io, ias) &
				  +conjg(zil(l2)) * riloa(ilo, l1, io, l2, ias, igq) * &
				  conjg(ylmgq(lm2, igq, iq))* &
				  xsgnt(lm1, lm2, lm3)
			  end do
		       end do
		       if (input%xs%dbglev.gt.2) then
			  write(u2, '(6i5, 2g18.10)') &
			       igq, ias, m3, ilo, lm1, io, &
			       intrgalo(m3, ilo, lm1, io, ias)
		       end if
		    end do
		 end do
	      end do
	   end do
	end do
        !-------------------------------------!
        !     local-orbital-APW integrals     !
        !-------------------------------------!
	do ilo=1, nlorb(is)
	   l1=lorbl(ilo, is)
	   do cm1=1, m1shape(l1)
	      m1=m1map(l1, cm1)
	      lm1=idxlm(l1, m1)
	      do cl2=1, l2shape(l1, m1)
		 l3=l2map(l1, m1, cl2)
		 do cm2=1, m2shape(l1, m1, l3)
		    m3=m2map(l1, m1, l3, cm2)
		    lm3=idxlm(l3, m3)
		    do io=1, apword(l3, is)
		       do cl3=1, l3shape(l1, m1, l3, m3)
			  l2=l3map(l1, m1, l3, m3, cl3)
			  do cm3=1, m3shape(l1, m1, l3, m3, l2)
			     m2=m3map(l1, m1, l3, m3, l2, cm3)
			     lm2=idxlm(l2, m2)
			     intrgloa(m1, ilo, lm3, io, ias)= &
				  intrgloa(m1, ilo, lm3, io, ias) &
				  +conjg(zil(l2)) * riloa(ilo, l3, io, l2, ias, igq) * &
				  conjg(ylmgq(lm2, igq, iq))* &
				  xsgnt(lm1, lm2, lm3)
			  end do
		       end do
		       if (input%xs%dbglev.gt.2) then
			  write(u3, '(6i5, 2g18.10)') &
			       igq, ias, m1, ilo, lm3, io, &
			       intrgloa(m1, ilo, lm3, io, ias)
		       end if
		    end do
		 end do
	      end do
	   end do
	end do
        !-----------------------------------------------!
        !     local-orbital-local-orbital integrals     !
        !-----------------------------------------------!
	do ilo1=1, nlorb(is)
	   l1=lorbl(ilo1, is)
	   do cm1=1, m1shape(l1)
	      m1=m1map(l1, cm1)
	      lm1=idxlm(l1, m1)
	      do ilo2=1, nlorb(is)
		 l3=lorbl(ilo2, is)
		 do cm2=1, m2shape(l1, m1, l3)
		    m3=m2map(l1, m1, l3, cm2)
		    lm3=idxlm(l3, m3)
		    do cl3=1, l3shape(l1, m1, l3, m3)
		       l2=l3map(l1, m1, l3, m3, cl3)
		       do cm3=1, m3shape(l1, m1, l3, m3, l2)
			  m2=m3map(l1, m1, l3, m3, l2, cm3)
			  lm2=idxlm(l2, m2)
			  intrglolo(m1, ilo1, m3, ilo2, ias)= &
			       intrglolo(m1, ilo1, m3, ilo2, ias) &
			       +conjg(zil(l2)) * rilolo(ilo1, ilo2, l2, ias, igq)* &
			       conjg(ylmgq(lm2, igq, iq))* &
			       xsgnt(lm1, lm2, lm3)
		       end do
		    end do
		    if (input%xs%dbglev.gt.2) then
		       write(u4, '(6i5, 2g18.10)') &
			    igq, ias, m1, ilo1, m3, ilo2, &
			    intrglolo(m1, ilo1, m3, ilo2, ias)
		    end if
		 end do
	      end do
	   end do
	end do
        ! end loops over atoms and species
     end do
  end do
  ! deallocate
  if (input%xs%dbglev.gt.2) then
     ! close files
     close(u1)
     close(u2)
     close(u3)
     close(u4)
  end if
end subroutine ematgntsum
