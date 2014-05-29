# netJina
Tools for working with [JINA reaclib](https://groups.nscl.msu.edu/jina/reaclib/db/) and [starlib](http://starlib.physics.unc.edu).

## Dependencies

*   The code is written in FORTRAN, and uses the `utils` and `const` modules from [MESA](http://mesa.sourceforge.net). Compilation and testing have been done with the MESA SDK and with a patched version v6208 of MESA.

*   The code requires three databases, containing the nuclide database, the reaclib database, and the starlib database.  The script `fetch_data` fetches compresed databases from https://dl.dropboxusercontent.com/u/52649885/netJina/.

## Installation

### Prerequisie
Unfortunately, you need to patch the `utils` module of v6208 of MESA. This necessity will disappear once the latest version of MESA is released.

1.  Copy `utils_patch` into the `utils` directory of your MESA v6022 tree.
2.  Go to that directory: `cd $MESA_DIR/utils`
2.  Run `patch < utils_patch`
3.  Remake the library: `./i1`

### Basic

    ./install [-d /path/to/project/root]

This top-level script calls all of the other build scripts, and installs them in subdirectories of `/path/to/project/root`. The fortran `*.mod` include files go into a subdirectory `include`, while the compiled library goes into a subdirectory `lib`. If `install` is called with no options, then the files are installed in the local directory in `include` and `lib`.  The database files and any cache files are placed in a subdirectory `data`.

### Detailed

For a more granular build, first `./fetch_data`: if the data files are missing or have the incorrect checksum, then fresh copies are downloaded and checked. Then `./build_and_test`: the libraries are compiled, and a small test program is compiled and run.  The output of the test program is in `test/test_output` and should be compared to `test/sample_output`. The `*.a` and `*.mod` files should be copied from `make` into their final locations, as should `starlib_db`, `reaclib_db`, and `nuclib_db` from the `data` directory. The `data/cache` directory and contents should also be copied to its final location.

## How to use (still under development)
Look in `test/src/test_io.f` for an example of source code, and in `test/make/makefile` for an example of compiling and linking the libraries.

In summary, you initialize the module

    call netJina_init(datadir,nuclib,nuclide_dict,reaclib,starlib,rates_dict,starlib_dict,ierr)
    if (ierr /= 0) then
        write(error_unit,*) 'failure in initialization ',ierr
        stop
    end if

where
    
    character(len=*), parameter :: datadir = '../data'
    type(reaclib_data) :: reaclib
    type(starlib_data) :: starlib
    type(nuclib_data) :: nuclib
    integer :: ierr
    type(integer_dict), pointer :: rates_dict=>null(), nuclide_dict=>null()
    
Then, to get the reaction parameters

    write (output_unit,'(/,a)') 'What are the reaction channels for ca37?'
    iso = 'ca37'
    call make_channel_handles(iso,nuclib,nuclide_dict,handles,ierr)

This returns the handles -- identifiers for the reactions -- for the following channels.

    character(len=max_id_length),dimension(N_bdat_channels) :: handles

	i_pg							:= (p,g)
	i_an							:= (a,n)
	i_ag							:= (a,g)
	i_ap							:= (a,p)
	i_ng							:= (n,g)
	i_np							:= (n,p)
	i_gp							:= (g,p)
    i_na                            := (n,a)
	i_ga							:= (g,a)
	i_pa							:= (p,a)
	i_gn							:= (g,n)
	i_pn							:= (p,n)

Now that you have the handles, you can get the reaction parameters.  From reaclib:

    write(output_unit,'(/,a)') 'What are the reaction parameters for ca37'
    call get_bdat_channels(reaclib,rates_dict,handles,n_coeff,rate_coefficients,q,rate_mask)

    integer, dimension(N_bdat_channels) :: n_coeff
    real(dp), dimension(ncoefficients*max_terms_per_rate,N_bdat_channels) :: rate_coefficients
    real(dp), dimension(N_bdat_channels) :: q
    logical, dimension(N_bdat_channels) :: rate_mask

From starlib:

    write(output_unit,'(/,a)') 'with starlib...'
    call get_bdat_rates(starlib,starlib_dict,handles,T9,rate,uncertainty,q,rate_mask)

    real(dp), dimension(number_starlib_temps,N_bdat_channels) :: T9, rate, uncer
    tainty
    real(dp), dimension(N_bdat_channels) :: q
    logical, dimension(N_bdat_channels) :: rate_mask

Becauase reaclib has a standard seven-coefficient fit per rate, some reactions have multiple rate entries: the total reaction rate is the sum of the individual terms.  Thus, for the channel with index `channel_index`, `rate_coefficients(1:7,channel_index)` has the coefficients for the first entry, `rate_coefficients(8:14,channel_index)` has the second, with a total of `7*n_coeff(channel_index)` coefficients for that reactions.  The Q-values are stored in `q(channel_index)`. Finally, logical array `rate_mask` contains for each channel the value `.TRUE.` if the rate is present in reaclib and is a forward rate.
	
## How it works
After reading in the `reaclib` database, the code generates for each rate a "handle" and uses that to build a dictionary.  When called with an (isotope, channel) pair a handle is easily constructed and used to look up the HEAD rate in the database. The number of entries for that reaction is stored in the database, so the code just needs to read those terms and pack them into the output.

## To do
*   The generation of handles and returning of reaction parameters should be handled by one wrapper routine, that would also take care of exceptions.
