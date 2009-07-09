



subroutine readbroydsteps_and_init_SY(noldsteps, n, S, Y, potential, residual)

	use modmixermsec, only: record_of_last_iter, noldstepsmax
	implicit none
	integer, intent(in)::noldsteps, n
	real(8), INTENT(OUT)::S(n, noldstepsmax), Y(n, noldstepsmax)
	real(8), INTENT(in)::potential(n), residual(n)
	integer::i, skipp
	character(256), external:: outfilenamestring
	character(256)::filetag
	integer::reclength, rectoread, firstrec
	inquire(iolength=reclength) potential, residual
	filetag="BROYDEN"
	open(23, file = outfilenamestring(filetag, 1), ACCESS = "DIRECT", RECL = reclength, &
			ACTION= "READ", FORM = 'UNFORMATTED')
	if(noldsteps .lt. noldstepsmax) then
		firstrec=1
	else
		firstrec=mod(record_of_last_iter, noldstepsmax)+1
	endif
	S=0
	Y=0
	skipp =noldstepsmax-noldsteps
	Do i=1, noldsteps
		rectoread=firstrec-1+i
		if (rectoread.gt.noldstepsmax) rectoread=rectoread-noldstepsmax
		read(23, rec=rectoread) s(:, i+skipp), y(:, i+skipp)
	end do
	close(23)

	Do i=1, noldsteps
		S(:, i+skipp)=S(:, i+skipp)-potential
		Y(:, i+skipp)=Y(:, i+skipp)-residual
	end do


end subroutine
