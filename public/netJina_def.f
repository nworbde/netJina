!   netJina_def.f
!
!   defines data layout for the JINA reaclib database
!   Edward Brown, Michigan State University
!
!   Requires installation of MESA (mesa.sourceforge.net) utils lib, v 6022
!

module netJina_def
    use utils_def, only: integer_dict

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
		integer,dimension(:),allocatable :: chapter
		character(len=iso_name_length),dimension(:,:),allocatable :: species
		character(len=iso_name_length),dimension(:),allocatable :: label
		character,dimension(:),allocatable :: reaction_flag
		character,dimension(:),allocatable :: reverse_flag
		real,dimension(:),allocatable :: Qvalue
		real,dimension(:,:),allocatable :: coefficients
        type(integer_dict), pointer :: reaclib_dict
	end type reaclib_data

	! container to hold locations for all terms for a given rate. Reaclib uses a
	! seven-term fit for each entry.  If this is insufficient, it uses more &
	! than one entry per reaction.
	type rate_location
		character(len=max_id_length) :: reaction_handle
		integer :: nterms
		integer, dimension(max_terms_per_rate) :: indices
	end type rate_location

end module netJina_def
