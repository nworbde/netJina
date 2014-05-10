module netJina_bdat

contains
    
    subroutine do_get_bdat_channels(reaclib,rates_dict, &
    & handles,n_coeff,rate_coefficients,q,rate_mask)
        use utils_def, only: integer_dict
        use utils_lib, only: integer_dict_lookup
        use netJina_def
        
        type(reaclib_data), intent(in) :: reaclib
        type(integer_dict), pointer :: rates_dict
        character(len=max_id_length), intent(in),  &
        & dimension(N_bdat_channels) :: handles
        integer, intent(out), dimension(N_bdat_channels) :: n_coeff
        real(dp), dimension(ncoefficients*max_terms_per_rate,N_bdat_channels), &
        & intent(out) :: rate_coefficients
        real(dp), intent(out), dimension(N_bdat_channels) :: q
        logical, intent(out), dimension(N_bdat_channels) :: rate_mask
        integer :: ierr
        integer :: i, indx, istart, iend

        ierr = 0
        do i = 1, N_bdat_channels
            call integer_dict_lookup(rates_dict,trim(handles(i)),indx,ierr)
            if (ierr /= 0) then
                rate_mask(i) = .FALSE.
                cycle
            end if
            q(i) = reaclib% Qvalue(indx)
            n_coeff(i) = reaclib% N_rate_terms(indx)*ncoefficients
            istart = indx
            iend = istart + reaclib% N_rate_terms(indx) - 1
            rate_coefficients(1:n_coeff(i),i) =  reshape( &
            & reaclib% coefficients(1:ncoefficients,istart:iend),[n_coeff(i)])
            rate_mask(i) = (reaclib% reverse_flag(indx) /= 'v')
        end do
    end subroutine do_get_bdat_channels
    
    subroutine do_make_channel_handles(isotope,nuclib,nuclide_dict,handles,ierr)
        use netJina_def
        use utils_def
        use utils_lib
        use iso_fortran_env, only: error_unit
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
        character(len=*), parameter :: null_handle = 'None'
        
        ! lookup nuclide
        call integer_dict_lookup(nuclide_dict,trim(isotope),indx,ierr)
        if (failure('looking for '//trim(isotope)//'. error code ',ierr)) &
        & return
        Z = nuclib% Z(indx)
        N = nuclib% N(indx)

        handles = null_handle
        
        ! set limits on Z, N
        if (Z <= 2 .or. Z > max_element_z-2) then
            ierr = -10
            write(error_unit,'(a,i3,a)') 'Z = ',Z, &
            & ' will produce out-range-links'
            return
        end if
        if (N <= 2 .or. N > maxval(nuclib% N)-2) then
            ierr = -11
            write(error_unit,'(a,i3,a)') 'N = ',N, &
            & ' will produce out-of-range links'
            return
        end if
        
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
