!   netJina_def.f
!
!   defines data layout for the JINA reaclib database
!   Edward Brown, Michigan State University
!
!   Requires installation of MESA (mesa.sourceforge.net) v 6022 for utils and 
!   const modules
!

module netJina_def
    use const_def, only: sp,dp

    ! data storage parameter for nuclib
    ! each isotope has a name and "provenance" -- a reference to where the 
    ! mass, spin, partition function are described.
    integer, parameter :: iso_name_length = 5
    integer, parameter :: provenance_length = 6
    integer, parameter :: reaction_reference_length = 4
    ! maximum number of nuclides in nucchem database
    integer, parameter :: max_nnuclib=10000
    ! no. entries in partition fcn table
    integer, parameter :: npfcn = 24

    ! maximum number of individual rates in reaclib.
    integer, parameter :: max_nreaclib=120000    
    ! data storage parameters for reaclib
    integer, parameter :: max_species_per_reaction=6
    integer, parameter :: ncoefficients=7
    integer, parameter :: nchapters=11
    integer, parameter :: max_terms_per_rate = 12
    ! we keep track of reactions using a unique "handle" generated from the 
    ! isotope names
    integer, parameter :: max_id_length = 42
    ! max. Z in the nucchem database
    integer, parameter :: max_element_z = 112
    ! used to output a formatted reaction string
    integer, parameter :: length_reaction_string = 45

    ! for starlib
    ! maximum number of individual rates in starlib
    integer, parameter :: max_nstarlib = 80000
    ! data storage parameters for starlib
    integer, parameter :: number_starlib_temps = 60
    integer, parameter :: starlib_max_species_per_reaction = 6

    
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
    
    ! storage for nuclib database
    type nuclib_data
        integer :: Nnuclides
        ! for each nuclide, store its name, provenance, mass (real), charge and 
        ! neutron numbers, ground-state spin, mass excess (MeV) and partition 
        ! function table
        character(len=iso_name_length), dimension(:), allocatable :: name
        character(len=provenance_length),dimension(:), allocatable :: provenance
        real(dp), dimension(:), allocatable :: A
        integer, dimension(:), allocatable  :: Z
        integer, dimension(:), allocatable :: N
        real(dp), dimension(:), allocatable :: spin
        real(dp), dimension(:), allocatable :: mass_excess
        real(dp), dimension(:,:), allocatable :: pfcn
    end type nuclib_data
    
    ! temperatures, partition function evaluation
    real(dp), dimension(npfcn) :: pfcn_T9
    
    ! flags for reaction channels, bdat-style
    integer, parameter :: i_pg = 1
    integer, parameter :: i_an = i_pg+1
    integer, parameter :: i_ag = i_an+1
    integer, parameter :: i_ap = i_ag+1
    integer, parameter :: i_ng = i_ap+1
    integer, parameter :: i_np = i_ng+1
    integer, parameter :: i_gp = i_np+1
    integer, parameter :: i_na = i_gp+1
    integer, parameter :: i_ga = i_na+1
    integer, parameter :: i_pa = i_ga+1
    integer, parameter :: i_gn = i_pa+1
    integer, parameter :: i_pn = i_gn+1
    integer, parameter :: N_bdat_channels = i_pn
    
    integer,dimension(N_bdat_channels) ::  &
    & bdat_dZ = [1,2,2,1,0,-1,-1,-2,-2,-1,0,1],  &
    & bdat_dN = [0,1,2,2,1,1,0,-1,-2,-2,-1,-1]
    
    ! flags for reaction types: r_<number in>_<number out> = chapter id
    integer, parameter :: &
    &   r_one_one   = 1, &
    &   r_one_two   = 2, &
    &   r_one_three = 3, &
    &   r_two_one   = 4, &
    &   r_two_two   = 5, &
    &   r_two_three = 6, &
    &   r_two_four  = 7, &
    &   r_three_one = 8, &
    &   r_three_two = 9, &
    &   r_four_two  = 10,&
    &   r_one_four  = 11
    
    ! Nin(chapter), Nout(chapter) give number of nuclides on the entrance and 
    ! exit channels for that chapter
    integer, dimension(nchapters) :: nJ_Nin = (/1,1,1,2,2,2,2,3,3,4,1/)
    integer, dimension(nchapters) :: nJ_Nout = (/1,2,3,1,2,3,4,1,2,2,4/)
    
    ! storage containers for rate data
    type rate_data
        integer :: Nentries
        integer,dimension(:),allocatable :: chapter
        character(len=iso_name_length),dimension(:,:),allocatable :: species
        character(len=reaction_reference_length),dimension(:),allocatable :: label
        character,dimension(:),allocatable :: reaction_flag
        character,dimension(:),allocatable :: reverse_flag
        real(dp),dimension(:),allocatable :: Qvalue
    end type rate_data
    
    ! storage container for reaclib file
    type, extends(rate_data) :: reaclib_data
        real(dp),dimension(:,:),allocatable :: coefficients
        integer,dimension(:),allocatable :: N_rate_terms
    end type reaclib_data
    
    ! storage container for starlib file
    type, extends(rate_data) :: starlib_data
        real(dp),dimension(:,:), allocatable :: T9
        real(dp),dimension(:,:), allocatable :: rate
        real(dp),dimension(:,:), allocatable :: uncertainty
    end type starlib_data
    
end module netJina_def
