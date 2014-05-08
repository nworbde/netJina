program test_io
    use iso_fortran_env, only: output_unit, error_unit
    use netJina_def
    use netJina_lib
    use utils_lib
    
    character(len=*), parameter :: reaclib_file = '../data/20140507default2'
    type(reaclib_data) :: reaclib
    integer :: ierr
    integer :: i,j
    type(integer_dict), pointer :: rates_dict=>null()
    integer :: max_terms, indx(1), dindx
    character(len=max_id_length) :: handle
    
    call load_reaclib(reaclib_file,reaclib,rates_dict,ierr)
    
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
    max_terms = maxval(reaclib% N_rate_terms)
    indx = maxloc(reaclib% N_rate_terms)
    
    write(output_unit,*) 'reaction: ',reaclib% species(:,indx(1))
    write(output_unit,*) 'has ',max_terms,' terms'
    
    call get_handle(reaclib,indx(1),handle)
    call integer_dict_lookup(rates_dict,handle,dindx,ierr)
    if (ierr == 0) then
        print *, 'rates_dict returns ',dindx,' for rate head.'
        do i = dindx, dindx+max_terms-1
            print *,reaclib% species(:,i),reaclib% N_rate_terms(i)
        end do
    end if
end program test_io
