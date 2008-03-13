module jacobidavidsoncommon
use modfvsystem
type(evsystem)::system
complex(8),allocatable ::p(:)
integer,allocatable::ipiv(:)

end module