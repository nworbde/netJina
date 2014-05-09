module netJina_bdat

contains
    subroutine do_make_channel_handles(isotope,nuclib,nuclide_dict,handles,ierr)
        use netJina_def
        use utils_def
        use utils_lib
        character(len=iso_name_length), intent(in) :: isotope
        type(nuclib_data), intent(in) :: nuclib
        type(integer_dict), pointer :: nuclide_dict
        character(len=max_id_length),dimension(N_bdat_channels),intent(out) :: &
        & handles
        integer, intent(out) :: ierr
        character(len=2),parameter :: nj = 'nJ'
        character(len=3),parameter :: sep = '_to'
        character(len=5),parameter,dimension(N_bdat_channels) ::  &
        & lhs = [character(len=5) :: '_p_','_he4_','_he4_','_he4_','_n_', &
        & '_n_','_','_n_','_','_p_','_','_p_']
        character(len=5),parameter,dimension(N_bdat_channels) ::  &
        & rhs = [character(len=5) :: '_','_n_','_','_p_','_','_p_','_p_', &
        & '_he4_','_he4_','_he4_','_n_','_n_']
        integer :: indx
        character(len=iso_name_length), dimension(N_bdat_channels) :: products
        integer :: Z, N, Zt, Nt, i
        
        ! lookup nuclide
        call integer_dict_lookup(nuclide_dict,trim(isotope),indx,ierr)
        if (failure('looking for '//trim(isotope)//'. error code ',ierr)) &
        & return
        Z = nuclib% Z(indx)
        N = nuclib% N(indx)
        
        do i = 1, N_bdat_channels
            Zt = Z + bdat_dZ(i)
            Nt = N + bdat_dN(i)
            write(products(i),'(a,i0)') trim(element_name(Zt)),Zt+Nt
        end do
        
        do i = 1, N_bdat_channels
            handles(i) = nj//trim(lhs(i))//trim(isotope)// &
            & sep//trim(rhs(i))//trim(products(i))
        end do
        
    end subroutine do_make_channel_handles

    function failure(action,ierr)
        use iso_fortran_env, only: error_unit
        character(len=*), intent(in) :: action
        integer, intent(in) :: ierr
        logical :: failure
        character(len=*), parameter :: err_format = &
        &  '("Error while ",a,". error code = ",i0)'
        
        if (ierr == 0) then
            failure = .FALSE.
            return
        end if
        failure = .TRUE.
        write (error_unit,err_format) action, ierr
    end function failure
    

end module netJina_bdat
