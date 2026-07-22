!> This preprocessor file controls the inclusion of assertion checks.
!> If USE_ASSERT is defined, it includes the 'asserts' module and defines
!> the macro CALL_ASSERT to expand to a call to the assert function. Otherwise,
!> CALL_ASSERT is replaced by a comment, effectively disabling assertions.
!>
!> This approach ensures that assert calls are completely omitted from
!> the compiled code when disabled, even if the condition involves a
!> non-pure function, because they are replaced by comments.

#if defined(USE_ASSERT)
use asserts, only: assert
#define CALL_ASSERT call assert
#else
#define CALL_ASSERT !!
#endif

