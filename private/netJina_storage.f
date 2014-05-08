module netJina_storage

contains
	subroutine allocate_reaclib_data(r, n, ierr)
        use netJina_def
		type(reaclib_data), intent(out) :: r
		integer, intent(in) :: n
		integer, intent(out) :: ierr
	
		ierr = 0
		allocate(r%chapter(n),r%species(max_species_per_reaction,n), &
		& r%label(n),r%reaction_flag(n), &
		& r%reverse_flag(n),r%Qvalue(n),r%coefficients(ncoefficients,n), &
		& r% N_rate_terms(n), stat=ierr)
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

end module netJina_storage