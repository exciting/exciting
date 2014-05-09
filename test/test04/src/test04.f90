program test
	use modmixermsec,only:residual,last_outputp,initmixermsec,&
	freearraysmixermsec,noldstepsmax,noldsteps
	use modreport
	use modmain,only:CHGIR
	implicit none
	real(8),allocatable::S(:,:),Y(:,:),YY(:,:),STEP(:),Potential(:)
	integer::iscl,n,teststeps,steps,i
	logical::passed1,passed2
	passed1=.true.
	passed2=.true.
	CHGIR=1.0
		n=4
		teststeps=24
        noldstepsmax=8
	allocate (S(n,noldstepsmax),Y(n,noldstepsmax),YY(noldstepsmax,noldstepsmax)&
	,STEP(n))
	allocate (Potential(n))

	call initmixermsec(n,noldstepsmax)
	do iscl=1,teststeps
		potential=4*iscl
		residual=2*iscl
		call write_current_to_broyden_file(n,iscl,potential,residual)
	end do
		call readbroydsteps_and_init_SY(noldsteps,n,S,Y,potential,residual)
	do iscl=1,noldsteps
		steps=teststeps-noldsteps+iscl
		write(*,*)S(:,iscl),"=",(4*steps-4*teststeps)
		write(*,*)Y(:,iscl),"=",(2*steps-2*teststeps)
		do i=1,n
			if(S(i,iscl).ne.(4*steps-4*teststeps))passed1=.false.
			if(Y(i,iscl).ne.(2*steps-2*teststeps))passed2=.false.
		end do
	end do
testplan_name ="test4"
tdirectory="test04"
	call inittestoutputfile(50) !file unit 50

	testunitname="Broydenfile S rwtest"
	inputf="-"
	outputf="-"
	call testreport(passed1)

	testunitname="Broydenfile Y rwtest"
	inputf="-"
	outputf="-"

	call testreport(passed2)

	call finalizeoutput(50)

end program
