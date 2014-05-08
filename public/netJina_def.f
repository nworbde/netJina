!   netJina_def.f
!
!   defines data layout for the JINA reaclib database
!   Edward Brown, Michigan State University
!
!   Requires installation of MESA (mesa.sourceforge.net) utils lib, v 6022
!

module netJina_def
    use const_def, only: sp,dp
    use utils_def, only: integer_dict
    use utils_lib, only: integer_dict_free

    ! reaclib uses a character handle for each isotope
    integer, parameter :: iso_name_length = 5
    ! maximum number of individual rates in reaclib.
	integer, parameter :: max_nreaclib=120000
	integer, parameter :: max_species_per_reaction=6
	integer, parameter :: ncoefficients=7
	integer, parameter :: nchapters=11
	integer, parameter :: ninverse_coeff = 2
	integer, parameter :: max_ec_rates = 8
	integer, parameter :: max_terms_per_rate = 20
    ! we keep track of reactions using a unique "handle" generated from the 
    ! isotopie names
	integer, parameter :: max_id_length = 42
    
	! flags for reaction types: r_<number in>_<number out> = chapter id
	integer, parameter :: &
	&	r_one_one	= 1, &
	&	r_one_two 	= 2, &
	&   r_one_three	= 3, &
	&   r_two_one	= 4, &
	&	r_two_two	= 5, &
	&   r_two_three	= 6, &
	&	r_two_four	= 7, &
	&   r_three_one	= 8, &
	&	r_three_two	= 9, &
	&   r_four_two	= 10,&
	&   r_one_four	= 11
    
    ! Nin(chapter), Nout(chapter) give number of nuclides on the entrance and 
    ! exit channels for that chapter
	integer, dimension(nchapters) :: Nin = (/1,1,1,2,2,2,2,3,3,4,1/)
	integer, dimension(nchapters) :: Nout = (/1,2,3,1,2,3,4,1,2,2,4/)
    
	! storage container for reaclib file
	type reaclib_data
        integer :: Nentries
		integer,dimension(:),allocatable :: chapter
		character(len=iso_name_length),dimension(:,:),allocatable :: species
		character(len=iso_name_length),dimension(:),allocatable :: label
		character,dimension(:),allocatable :: reaction_flag
		character,dimension(:),allocatable :: reverse_flag
		real(dp),dimension(:),allocatable :: Qvalue
		real(dp),dimension(:,:),allocatable :: coefficients
        type(integer_dict), pointer :: reaclib_dict
	end type reaclib_data

	! container to hold locations for all terms for a given rate. Reaclib uses a
	! seven-term fit for each entry.  If this is insufficient, it uses more
	! than one entry per reaction.
	type rate_location
		character(len=max_id_length) :: reaction_handle
		integer :: nterms
		integer, dimension(max_terms_per_rate) :: indices
	end type rate_location

contains
	subroutine allocate_reaclib_data(r, n, ierr)
        use utils_lib, only: integer_dict_free
		type(reaclib_data), intent(out) :: r
		integer, intent(in) :: n
		integer, intent(out) :: ierr
	
		ierr = 0
		allocate(r%chapter(n),r%species(max_species_per_reaction,n), &
		& r%label(n),r%reaction_flag(n), &
		& r%reverse_flag(n),r%Qvalue(n),r%coefficients(ncoefficients,n), &
		& stat=ierr)
!         if (associated(r% reaclib_dict)) then
!             call integer_dict_free(r% reaclib_dict)
            nullify(r% reaclib_dict)
!         end if
        r% Nentries = n
	end subroutine allocate_reaclib_data

	subroutine free_reaclib_data(r)
		type(reaclib_data), intent(inout) :: r
		if (allocated(r% chapter)) deallocate(r% chapter)
		if (allocated(r% species)) deallocate(r% species)
		if (allocated(r% label)) deallocate(r% label)
		if (allocated(r% reaction_flag)) deallocate(r% reaction_flag)
		if (allocated(r% reverse_flag)) deallocate(r% reverse_flag)
		if (allocated(r% Qvalue)) deallocate(r% Qvalue)
		if (allocated(r% coefficients)) deallocate(r% coefficients)
        if (associated(r% reaclib_dict)) call integer_dict_free(r% reaclib_dict)
        r% Nentries = 0
	end subroutine free_reaclib_data
    
    subroutine copy_reaclib_data(old_data,new_data,ierr)
        type(reaclib_data), intent(in) :: old_data
        type(reaclib_data), intent(out) :: new_data
        integer, intent(out) :: ierr
        integer :: n
        
        if (new_data% Nentries /= 0) call free_reaclib_data(new_data)
        call allocate_reaclib_data(new_data,old_data% Nentries, ierr)
        if (ierr /= 0) return
        
        n = old_data% Nentries
        new_data% chapter(:) = old_data% chapter(1:n)
		new_data% species(:,:) = old_data% species(:,1:n)
		new_data% label(:) = old_data% label(1:n)
		new_data% reaction_flag(:) = old_data% reaction_flag(1:n)
		new_data% reverse_flag(:) = old_data% reverse_flag(1:n)
		new_data% Qvalue(:) = old_data% Qvalue(1:n)
		new_data% coefficients(:,:) = old_data% coefficients(:,1:n)
    end subroutine copy_reaclib_data
    
end module netJina_def
