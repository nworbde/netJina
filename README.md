# netJina
Tools for working with [JINA reaclib](https://groups.nscl.msu.edu/jina/reaclib/db/).

## Dependencies

*   The code will be a FORTRAN module, and will make use of the `utils` module from [MESA](http://mesa.sourceforge.net).

*   The code requires two databases, one containing the nuclide database and the other containing the reaclib database.  These are obtained from the JINA website.  The script `install_data` does this automatically using `curl`.

## Installation

    ./install_data
    ./build_and_test

## How to use (still under development)
To use, in the top-level code

    use netJina_def
    use netJina_lib
    
    type(reaclib_data) :: reaclib
    type(nuclib_data) :: nuclib
    integer :: ierr
    type(integer_dict), pointer :: rates_dict=>null(), nuclide_dict=>null()
    
    call netJina_init(datadir,nuclib,nuclide_dict,reaclib,rate_dict,ierr)
        
        character(len=*), intent(in) :: datadir
        ! directory holding the two database files
        
        type(nuclib_data), intent(out) :: nuclib
        ! data storage for nuclide database
        
        type(integer_dict), pointer :: nuclide_dict
        ! dictionary, used to get index of rate
        
        type(reaclib_data), intent(out) :: reaclib
        ! data storage for reaclib database
        
        type(integer_dict), pointer :: rate_dict
        ! dictionary: when called with a handle, it returns the index of the first entry for that rate.
        ! (reaclib may have more than one entry for a rate, in which case one sums over all rates
        ! to get the total.)
        
        integer, intent(out) :: ierr
        ! error flag. 0 on exit means success.

Then, to get a rate compatible with bdat-based networks

	call get_rate(nuclib,nuclide_dict,reaclib,rate_dict,isotope,rates,n_coeff,q,status,ierr)
	
where

	ratedb 							:= reaclib database read in during init
	isotope 						:= string containing the isotope, e.g., 'fe56' or 'na23'
	rates(max_coeff,max_channels) 	:= array of coeffiencts for rates
	n_coeff(max_channels)			:= number of coefficients for each channel
	q(max_channels)					:= q-value for each channel
	status(max_channels)			:= logical mask arrary.  stores FALSE if the channel is unavailable 
										or is not a forward rate
	ierr							:= error flag, 0 means all okay.
	
The channels (12 in all) are accessed via the following pointers

	i_pg							:= (p,g)
	i_ag							:= (a,g)
	i_ap							:= (a,p)
	i_ng							:= (n,g)
	i_np							:= (n,p)
	i_gp							:= (g,p)
	i_an							:= (a,n)
	i_ga							:= (g,a)
	i_pa							:= (p,a)
	i_gn							:= (g,n)
	i_pn							:= (p,n)
	
## How it works
After reading in the `reaclib` database, the code generates for each rate a "handle" and uses that to build a dictionary.  When called with an isotope, channel pair a handle is easily constructed and used to look up the HEAD rate in the database. The number of terms is stored in the database, so the code just needs to read those terms and pack them into the output.
