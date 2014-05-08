program test_bdat
    use iso_fortran_env, only: output_unit, error_unit
    use netJina_def
    use netJina_lib
    use utils_lib
    
    character(len=*), parameter :: datadir = '../data'
    type(reaclib_data) :: reaclib
    type(nuclib_data) :: nuclib
    integer :: ierr
    integer :: i,j
    type(integer_dict), pointer :: rates_dict=>null(), nuclide_dict=>null()
    character(len=max_id_length),dimension(N_bdat_channels) :: handles
    
    call netJina_init(datadir,nuclib,nuclide_dict,reaclib,rates_dict,ierr)
    if (ierr /= 0) then
        write(error_unit,*) 'load reaclib returned error ',ierr
        stop
    end if

    call make_channel_handles(' fe56',nuclib,nuclide_dict,handles,ierr)
    write (output_unit,'(i3,tr2,a)') (i,trim(handles(i)), i=1,N_bdat_channels)
end program test_bdat
