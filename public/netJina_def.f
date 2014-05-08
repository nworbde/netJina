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
	integer, parameter :: max_terms_per_rate = 20
    ! we keep track of reactions using a unique "handle" generated from the 
    ! isotope names
	integer, parameter :: max_id_length = 42
    ! max. Z in the nucchem database
	integer, parameter :: max_element_z = 112
    
    ! table of elements; note that hydrogen is referred to as 'p' in the table 
    ! and neutrons are referred to as 'n'
	character(len=iso_name_length), dimension(0:max_element_z) ::  &
	& element_name = [character(len=iso_name_length) :: & 
	& 'neut','h','he','li','be','b','c','n','o','f','ne', 'na','mg','al','si', &
	& 'p', 's', 'cl', 'ar', 'k', 'ca', 'sc', 'ti', 'v','cr','mn', 'fe', 'co', &
	& 'ni','cu','zn','ga','ge','as','se','br','kr','rb','sr','y','zr', &
	& 'nb','mo','tc','ru','rh','pd','ag','cd','in','sn','sb','te','i','xe', &
	& 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy','ho','er', &
	& 'tm','yb', 'lu','hf','ta','w','re','os','ir','pt','au','hg', &
	& 'tl','pb','bi','po','at','rn','fr','ra','ac','th','pa','u','np','pu', &
	& 'am','cm','bk','cf','es','fm','md','no','lr','rf','db','sg','bh','hs', &
	& 'mt','ds','rg','cn' ]
    
	! synonyms of hydrogen isotopes, index by mass number
	character(len=iso_name_length), dimension(3) ::  &
	& h_isotopes = [character(len=iso_name_length) :: 'p','d','t']

	! aluminum-26 isomers
	character(len=iso_name_length), dimension(2:3) ::  &
	& al_isomers = [character(len=iso_name_length) :: 'al-6','al*6']
    
    ! flags for reaction channels, bdat-style
	integer, parameter :: i_pg = 1
	integer, parameter :: i_ag = i_pg+1
	integer, parameter :: i_ap = i_ag+1
	integer, parameter :: i_ng = i_ap+1
	integer, parameter :: i_np = i_ng+1
	integer, parameter :: i_gp = i_np+1
	integer, parameter :: i_an = i_gp+1
	integer, parameter :: i_ga = i_an+1
	integer, parameter :: i_pa = i_ga+1
	integer, parameter :: i_gn = i_pa+1
	integer, parameter :: i_pn = i_gn+1
    integer, parameter :: N_bdat_channels = i_pn
    
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
	integer, dimension(nchapters) :: nJ_Nin = (/1,1,1,2,2,2,2,3,3,4,1/)
	integer, dimension(nchapters) :: nJ_Nout = (/1,2,3,1,2,3,4,1,2,2,4/)
    
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
        integer,dimension(:),allocatable :: N_rate_terms
	end type reaclib_data

	! container to hold locations for all terms for a given rate. Reaclib uses a
	! seven-term fit for each entry.  If this is insufficient, it uses more
	! than one entry per reaction.
	type rate_location
		character(len=max_id_length) :: reaction_handle
		integer :: nterms
		integer, dimension(max_terms_per_rate) :: indices
	end type rate_location
    
end module netJina_def
