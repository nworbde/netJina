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
    subroutine netJina_init(datadir,nuclib,nuclide_dict,reaclib,rate_dict,ierr)
        use iso_fortran_env, only: error_unit
        use netJina_def
        use netJina_io
        use netJina_storage
        use utils_lib, only: integer_dict_size
        character(len=*),parameter :: reaclib_db = 'reaclib_db'
        character(len=*),parameter :: nuclib_db = 'nuclide_db'
        
        character(len=*), intent(in) :: datadir
        type(nuclib_data), intent(out) :: nuclib
        type(integer_dict), pointer :: nuclide_dict
        type(reaclib_data), intent(out) :: reaclib
        type(integer_dict), pointer :: rate_dict
        integer, intent(out) :: ierr

        character(len=160) :: reaclib_filename, nuclib_filename
        
        nuclib_filename = trim(datadir)//'/'//nuclib_db
        reaclib_filename = trim(datadir)//'/'//reaclib_db
        if (reaclib% Nentries /= 0) call free_reaclib_data(reaclib)

        ierr = 0
        write(error_unit,'(a)')  &
        & 'loading nuclib from '//trim(nuclib_filename)
        call do_load_nuclib(nuclib_filename,nuclib,ierr)
        write(error_unit,'(/,a,i0,a)')  &
        & 'done. ',nuclib% Nnuclides, &
        & ' nuclides retrieved. now writing nuclide dictionary...'
        call do_parse_nuclides(nuclib,nuclide_dict,ierr)
        write(error_unit,'(a//)') 'done.'
        
        write(error_unit,'(a)')  &
        & 'loading reaclib from '//trim(reaclib_filename)
        call do_load_reaclib(reaclib_filename,reaclib,ierr)
        write(error_unit,'(a,i0,a)') 'done. ',reaclib% Nentries, &
        & ' entries retrieved. now writing reaction dictionary...'
        call do_parse_rates(reaclib,rate_dict,ierr)
        write(error_unit,'(a,i0,a)') 'done. ',integer_dict_size(rate_dict), &
        & ' unique rates found.'
    end subroutine netJina_init

    subroutine get_handle(reaclib,indx,handle)
        use netJina_def, only: reaclib_data, nJ_Nin, nJ_Nout, max_id_length
        use netJina_io, only: do_generate_handle
        
        type(reaclib_data), intent(in) :: reaclib
        integer, intent(in) :: indx
        character(len=max_id_length), intent(out) :: handle
        handle = do_generate_handle(reaclib% species(:,indx), &
        & nJ_Nin(reaclib% chapter(indx)), nJ_Nout(reaclib% chapter(indx)))
    end subroutine get_handle
    
    subroutine make_channel_handles(isotope,nuclib,nuclide_dict,handles,ierr)
        use netJina_def
        use utils_def
        use utils_lib
        use netJina_bdat
        character(len=iso_name_length), intent(in) :: isotope
        type(nuclib_data), intent(in) :: nuclib
        type(integer_dict), pointer :: nuclide_dict
        character(len=max_id_length),dimension(N_bdat_channels),intent(out) :: &
        & handles
        integer, intent(out) :: ierr
        call do_make_channel_handles(isotope,nuclib,nuclide_dict,handles,ierr)
    end subroutine make_channel_handles
end module netJina_lib
