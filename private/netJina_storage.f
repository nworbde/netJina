module netJina_storage

contains
    
    subroutine allocate_nuclib_data(w,n,ierr)
        use netJina_def
        type(nuclib_data), intent(out) :: w
        integer, intent(in) :: n
        integer, intent(out) :: ierr
        
        ierr = 0
        w% Nnuclides = 0
        allocate(w% name(n), w% provenance(n), w% A(n), w% Z(n), w% N(n),  &
        & w% spin(n), w% mass_excess(n), w% pfcn(npfcn,n), stat=ierr)
        if (ierr /= 0) return
        w% Nnuclides = n
    end subroutine allocate_nuclib_data
    
    subroutine free_nuclib_data(w)
        use netJina_def
        type(nuclib_data), intent(inout) :: w
        if (allocated(w% name)) deallocate(w% name)
        if (allocated(w% provenance)) deallocate(w% provenance)
        if (allocated(w% A)) deallocate(w% A)
        if (allocated(w% Z)) deallocate(w% Z)
        if (allocated(w% N)) deallocate(w% N)
        if (allocated(w% spin)) deallocate(w% spin)
        if (allocated(w% mass_excess)) deallocate(w% mass_excess)
        if (allocated(w% pfcn)) deallocate(w% pfcn)
        w% Nnuclides = 0
    end subroutine free_nuclib_data
    
    subroutine allocate_reaclib_data(r, n, ierr)
        use netJina_def
        type(reaclib_data), intent(out) :: r
        integer, intent(in) :: n
        integer, intent(out) :: ierr
    
        ierr = 0
        r% Nentries = n
        allocate(r%chapter(n),r%species(max_species_per_reaction,n), &
        & r%label(n),r%reaction_flag(n), &
        & r%reverse_flag(n),r%Qvalue(n),r%coefficients(ncoefficients,n), &
        & r% N_rate_terms(n), stat=ierr)
        if (ierr /= 0) return
        r% Nentries = n
    end subroutine allocate_reaclib_data

    subroutine free_reaclib_data(r)
        use netJina_def
        type(reaclib_data), intent(inout) :: r
        if (allocated(r% chapter)) deallocate(r% chapter)
        if (allocated(r% species)) deallocate(r% species)
        if (allocated(r% label)) deallocate(r% label)
        if (allocated(r% reaction_flag)) deallocate(r% reaction_flag)
        if (allocated(r% reverse_flag)) deallocate(r% reverse_flag)
        if (allocated(r% Qvalue)) deallocate(r% Qvalue)
        if (allocated(r% coefficients)) deallocate(r% coefficients)
        if (allocated(r% N_rate_terms)) deallocate(r% N_rate_terms)
        r% Nentries = 0
    end subroutine free_reaclib_data
    
    subroutine allocate_starlib_data(r, n, ierr)
        use netJina_def
        type(starlib_data), intent(out) :: r
        integer, intent(in) :: n
        integer, intent(out) :: ierr
        
        ierr = 0
        allocate(r% chapter(n),  &
        &   r% species(starlib_max_species_per_reaction,n), &
        &   r% label(n), &
        &   r% reverse_flag(n), &
        &   r% Qvalue(n), &
        &   r% T9(number_starlib_temps,n), &
        &   r% rate(number_starlib_temps,n), &
        &   r% uncertainty(number_starlib_temps,n), stat=ierr)
        r% Nentries = n
    end subroutine allocate_starlib_data

    subroutine free_starlib_data(r)
        use netJina_def
        type(starlib_data), intent(inout) :: r
        if (allocated(r% chapter)) deallocate(r% chapter)
        if (allocated(r% species)) deallocate(r% species)
        if (allocated(r% label)) deallocate(r% label)
        if (allocated(r% reverse_flag)) deallocate(r% reverse_flag)
        if (allocated(r% Qvalue)) deallocate(r% Qvalue)
        if (allocated(r% T9)) deallocate(r% T9)
        if (allocated(r% rate)) deallocate(r% rate)
        if (allocated(r% uncertainty)) deallocate(r% uncertainty)
        r% Nentries = 0
    end subroutine free_starlib_data

end module netJina_storage