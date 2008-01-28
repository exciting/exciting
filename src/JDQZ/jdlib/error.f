      subroutine error (m)
c
c     Coded by Diederik R. Fokkema
c
c     $Id: error.f,v 1.4 1995/07/26 09:26:26 fokkema Exp $
c
c     .. Parameters ..
c
      implicit none
      character m*(*)

ctex@ \begin{manpage}{ERROR} 
ctex@ \subtitle{Name}
ctex@    ERROR --- Type an error message and stop
ctex@
ctex@ \subtitle{Declaration}
ctex@    %declaration
ctex@ \subtitle{Parameters}
ctex@    \variable{m}
ctex@       character string. On entry m must contain the message string.
ctex@
ctex@ \subtitle{Description}
ctex@    This subroutine types an error message and stops.
ctex@
ctex@ \end{manpage}
ctex@ \begin{verbatim}
ctex@    % actual code
ctex@ \end{verbatim}
c
c     .. Local ..
c
c     None.
c
c     .. Called subroutines
c
c     None.
c
c     .. Executable statements
c
      print 10, m
 10   format (/,1x,'Error: ',a,/)
c
c     --- Stop
c
      stop
      end
