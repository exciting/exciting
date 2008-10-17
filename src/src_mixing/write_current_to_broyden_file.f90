subroutine write_current_to_broyden_file(n,iscl,potential,residual)
	use modmixermsec,only: record_of_last_iter,noldstepsmax
	implicit none
	integer ,intent(in)::n,iscl
	real(8),intent(in)::potential(n),residual(n)
	integer::reclength

	record_of_last_iter=mod(record_of_last_iter,noldstepsmax)+1
	inquire(iolength=reclength) potential,residual

	open(23,file="BROYDEN.OUT",ACCESS="DIRECT",RECL=reclength,FORM='UNFORMATTED')
	write(23,rec=record_of_last_iter)potential,residual
	close(23)
end subroutine
