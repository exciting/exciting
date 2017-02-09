Subroutine ematsumuo (iq,ik, igq, integral, xi)
	Use modmain
	Use modinput
	Use modxs
	Use m_getunit 
	Use modxas
	Implicit none
	Integer, Intent (In) :: iq, ik, igq
	Complex(8), Intent (In) :: integral(input%xs%lmaxemat+1,lmmaxapw,nxas,sto2-sta2+1)
	Complex(8), Intent (InOut):: xi(xasstop-xasstart+1, sto2-sta2+1, ngq(iq))
	! local variables
	Integer :: n1, n2, l2, lmax2, m2, lm2, l3, m3, lm3, ias,ia, is
	Complex(8) :: prefactor
	Real (8) :: vk (3), vq(3), vkq(3)
	! Setting xioo to zero
	xi(:,:,igq)=zzero
	is=input%xs%bse%xasspecies
	ia=input%xs%bse%xasatom
	ias=idxas(ia,is)
	lmax2=input%xs%lmaxemat
	
	Do n1=1,nxas
		Do n2=sta2,sto2
			Do l2=0,lmax2
				Do m2=-l2,l2
					lm2=idxlm(l2,m2)
					Do l3=0, input%groundstate%lmaxapw
						Do m3=-l3,l3
							lm3=idxlm(l3,m3)
							xi(n1,n2,igq)=xi(n1,n2,igq)+conjg (zil(l2))*integral(l2+1,lm3,n1,n2)*conjg &
									& (ylmgq(lm2, igq, iq)) * xsgntuo &
									& (n1, lm2, lm3)
						End Do
					End Do
				End Do
			End Do
		End Do
	End Do
	vk(:) = vkl (:, ik)
	vq(:) = vql(:,iq)
	vkq(:)=vk(:)+vq(:)
	prefactor=fourpi*conjg(sfacgq(igq, ias, iq))
	xi(:,:,igq)=xi(:,:,igq)*prefactor
End Subroutine ematsumuo
