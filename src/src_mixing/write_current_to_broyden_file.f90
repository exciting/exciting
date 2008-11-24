subroutine write_current_to_broyden_file(n,iscl,potential,residual)
	use modmixermsec,only: record_of_last_iter,noldstepsmax,noldstepsin_file
	implicit none
	integer ,intent(in)::n,iscl
	real(8),intent(in)::potential(n),residual(n)
	integer::reclength
	character(256), external:: outfilenamestring
	character(256)::filetag
	filetag="BROYDEN"
	record_of_last_iter=mod(record_of_last_iter,noldstepsmax)+1
	inquire(iolength=reclength) potential,residual

	open(23,file=outfilenamestring(filetag,1),ACCESS="DIRECT",RECL=reclength,FORM='UNFORMATTED')
	write(23,rec=record_of_last_iter)potential,residual
	close(23)
	noldstepsin_file=noldstepsin_file+1
	noldstepsin_file= min(noldstepsin_file,noldstepsmax)
end subroutine
