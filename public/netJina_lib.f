!   netJina_lib.f
!
!   interface routines for the JINA reaclib database
!   Edward Brown, Michigan State University
!
!   Requires installation of MESA (mesa.sourceforge.net) utils lib, v 6022
!

module netJina_lib
    
contains
    subroutine load_reaclib(filename,reaclib,ierr)
        use netJina_def
        use reaclib_io
        character(len=*), intent(in) :: filename
        type(reaclib_data), intent(out) :: reaclib
        integer, intent(out) :: ierr
        
        ierr = 0
        if (reaclib% Nentries /= 0) call free_reaclib_data(reaclib)
        print *,'loading reaclib'
        call do_load_reaclib(filename,reaclib,ierr)
        print *,'reaclib loaded'
    end subroutine load_reaclib

end module netJina_lib
