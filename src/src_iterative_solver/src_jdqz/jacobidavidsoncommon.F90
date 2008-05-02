module jacobidavidsoncommon
	use modfvsystem
	type(evsystem)::system
	type(HermiteanMatrix)::p
	complex(8),allocatable::pdiag(:)
end module
