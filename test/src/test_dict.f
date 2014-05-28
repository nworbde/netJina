program test_dict
    use, intrinsic :: iso_fortran_env, only: error_unit, iostat_end
    use utils_def, only: integer_dict
    use utils_lib, only: integer_dict_define_and_check, integer_dict_lookup, integer_dict_free, integer_dict_create_hash
    
    type(integer_dict), pointer :: rate_dict=>null()
    integer :: i, ikey, indx, ierr, ios
    character(len=42) :: handle
    integer :: handle_unit
    logical :: duplicate
    integer :: n_rates
    
    open(newunit=handle_unit,file='handles',status='old',action='read')
    i = 1
    n_rates = 0
    do
        read(handle_unit,*,iostat=ios) handle
        if (ios == iostat_end) exit
        call integer_dict_define_and_check(rate_dict,trim(handle),i,duplicate,ierr)
        if (ierr /= 0) stop 'whoa!'
        if ( .not. duplicate) n_rates = n_rates + 1
        if (mod(i,1000)==0) write (error_unit,'(a)',advance='no') '.'
        i = i+1
    end do

    print *, i-1, ' rates parsed. ',n_rates,' non-duplicate rates found'
end program test_dict
