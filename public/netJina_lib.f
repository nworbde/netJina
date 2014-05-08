!   netJina_lib.f
!
!   interface routines for the JINA reaclib database
!   Edward Brown, Michigan State University
!
!   Requires installation of MESA (mesa.sourceforge.net) utils lib, v 6022
!

module netJina_lib
    use utils_def, only: integer_dict
    
contains
    subroutine load_reaclib(filename,reaclib,rate_dict,ierr)
        use iso_fortran_env, only: error_unit
        use netJina_def
        use reaclib_io
        character(len=*), intent(in) :: filename
        type(reaclib_data), intent(out) :: reaclib
        type(integer_dict), pointer :: rate_dict
        integer, intent(out) :: ierr
        
        ierr = 0
        if (reaclib% Nentries /= 0) call free_reaclib_data(reaclib)
        
        write(error_unit,*) 'loading reaclib from '//trim(filename)//'...'
        call do_load_reaclib(filename,reaclib,ierr)
        write(error_unit,*) 'done. writing reaction dictionary...'
        call do_parse_rates(reaclib,rate_dict,ierr)
        
    end subroutine load_reaclib

end module netJina_lib
