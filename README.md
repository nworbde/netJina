# netJina
Tools for working with [JINA reaclib](https://groups.nscl.msu.edu/jina/reaclib/db/).

## Calling sequence
The code will be a FORTRAN module, and will make use of the `utils` module from [MESA](http://mesa.sourceforge.net).

To use, in the top-level code

	call reaclib_init(filename,ratedb,ierr)

where

	! filename 						:= name of reaclib file
	! ierr 							:= error flag, 0 means init went okay.
	! ratedb 						:= datastructure holding the reaclib rates
	
to load the reaclib database.

Then, to get a rate compatible with bdat-based networks

	call get_rate(ratedb,isotope,rates,n_coeff,q,status,ierr)
	
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
