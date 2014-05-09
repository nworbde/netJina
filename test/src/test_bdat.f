program test_bdat
    use iso_fortran_env, only: output_unit, error_unit
    use netJina_def
    use netJina_lib
    use utils_lib
    
    character(len=*), parameter :: datadir = '../data'
    type(reaclib_data) :: reaclib
    type(nuclib_data) :: nuclib
    integer :: ierr
    integer :: i,j,i_rate
    type(integer_dict), pointer :: rates_dict=>null(), nuclide_dict=>null()
    character(len=max_id_length),dimension(N_bdat_channels) :: handles
    character(len=iso_name_length) :: iso
    
    call netJina_init(datadir,nuclib,nuclide_dict,reaclib,rates_dict,ierr)
    if (ierr /= 0) then
        write(error_unit,*) 'load reaclib returned error ',ierr
        stop
    end if
    iso = 'fe56'
    call make_channel_handles(iso,nuclib,nuclide_dict,handles,ierr)
    do i = 1, N_bdat_channels
        call integer_dict_lookup(rates_dict,trim(handles(i)),i_rate,ierr)
        if (ierr == 0) then
            write (output_unit,*) i,trim(handles(i)), &
            & (trim(reaclib% species(j,i_rate)),j=1,4)
        else
            write (output_unit,*) i,trim(handles(i)), 'rate not found'
        end if
    end do

end program test_bdat
