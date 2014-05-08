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
        use utils_lib, only: integer_dict_size
        
        character(len=*), intent(in) :: filename
        type(reaclib_data), intent(out) :: reaclib
        type(integer_dict), pointer :: rate_dict
        integer, intent(out) :: ierr
        
        ierr = 0
        if (reaclib% Nentries /= 0) call free_reaclib_data(reaclib)
        
        write(error_unit,'(a)') 'loading reaclib from '//trim(filename)//'...'
        call do_load_reaclib(filename,reaclib,ierr)
        write(error_unit,'(a,i0,a)') 'done. ',reaclib% Nentries, &
        & ' entries retrieved. now writing reaction dictionary...'
        call do_parse_rates(reaclib,rate_dict,ierr)
        write(error_unit,'(a,i0,a)') 'done. ',integer_dict_size(rate_dict), &
        & ' unique rates found.'
    end subroutine load_reaclib

    subroutine get_handle(reaclib,indx,handle)
        use netJina_def, only: reaclib_data, nJ_Nin, nJ_Nout, max_id_length
        use reaclib_io, only: generate_handle
        
        type(reaclib_data), intent(in) :: reaclib
        integer, intent(in) :: indx
        character(len=max_id_length), intent(out) :: handle
        handle = generate_handle(reaclib% species(:,indx), &
        & nJ_Nin(reaclib% chapter(indx)), nJ_Nout(reaclib% chapter(indx)))
    end subroutine get_handle
end module netJina_lib
