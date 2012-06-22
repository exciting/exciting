       subroutine setwpa(nomeg,npar,omega,iwpa)
!
! This subroutien set the freuency points used for the Pade's approximation in the case when 
! nomeg > npar
!
       implicit none 
       integer(4),intent(in):: nomeg,npar
       real(8),intent(in) :: omega
       integer(4),intent(out)::iwpa(npar)
        
       integer(4):: step,ipar

       step=nomeg/(npar-1)
       do ipar=1,npar-1
         iwpa(ipar)=(ipar-1)*step+1
       enddo
       iwpa(ipar)=nomeg

       endsubroutine 

