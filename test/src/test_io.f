program test_io
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
    integer :: max_terms, indx(1), dindx
    character(len=max_id_length) :: handle
    
    call netJina_init(datadir,nuclib,nuclide_dict,reaclib,rates_dict,ierr)
    if (ierr /= 0) then
        write(error_unit,*) 'load reaclib returned error ',ierr
        stop
    end if
    
    call integer_dict_lookup(nuclide_dict,'cn337',dindx,ierr)
    if (ierr /=0) then
        write(error_unit,*) 'cn337 is not on the list'
    else
        write(output_unit,*) nuclib% name(dindx), nuclib% A(dindx)
    end if
    
    do i = 1,reaclib% Nentries
        if (reaclib% chapter(i) == 5) exit
    end do
    if (reaclib% chapter(i) == 5) then
        write(output_unit,'(a)') 'Here are the first 50 entries of chapter 5'
        do j = i,i+50
            write(output_unit,'(a5,a5," ==> ",a5,a5)')  &
            & reaclib% species(1:4,j)
        end do
    else
        write(error_unit,'(a)') 'unable to find reactions in chapter 5'
    end if

    max_terms = maxval(reaclib% N_rate_terms)
    indx = maxloc(reaclib% N_rate_terms)
    
    write(output_unit,*) 'reaction: ',reaclib% species(:,indx(1)),  &
    & 'has ',max_terms,' terms'
    
    call get_handle(reaclib,indx(1),handle)
    call integer_dict_lookup(rates_dict,handle,dindx,ierr)
    if (ierr == 0) then
        do i = dindx, dindx+max_terms-1
            write(output_unit,'(i5,tr1,6a5,tr2,7es12.4)')  &
            & i,reaclib% species(:,i),reaclib% coefficients(:,i)
        end do
    end if
end program test_io
