program test_io
    use iso_fortran_env, only: output_unit, error_unit
    use netJina_def
    use netJina_lib
    
    character(len=*), parameter :: reaclib_file = '../data/20140507default2'
    type(reaclib_data) :: reaclib
    integer :: ierr
    integer :: i,j
    
    call load_reaclib(reaclib_file,reaclib,ierr)
    
    ! write the first 50 entries of chapter 5
    if (ierr == 0) then
        do i = 1,reaclib% Nentries
            if (reaclib% chapter(i) == 5) exit
        end do
        if (reaclib% chapter(i) == 5) then
            do j = i,i+50
                write(output_unit,'(a,a,"==>",a,a)') reaclib% species(1:4,j)
            end do
        else
            write(error_unit,'(a)') 'unable to find reactions in chapter 5'
        end if
    else
        write(error_unit,*) 'load reaclib returned error ',ierr
    end if
end program test_io
