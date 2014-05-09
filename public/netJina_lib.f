!   netJina_lib.f
!
!   interface routines for the JINA reaclib database
!   Edward Brown, Michigan State University
!
!   Requires installation of MESA (mesa.sourceforge.net) utils lib, v 6022
!

module netJina_lib
    use utils_def, only: integer_dict
    
contains
    subroutine netJina_init(datadir,nuclib,nuclide_dict,reaclib,rate_dict,ierr)
        use iso_fortran_env, only: error_unit
        use netJina_def
        use netJina_io
        use netJina_storage
        use utils_lib, only: integer_dict_size
        character(len=*),parameter :: reaclib_db = 'reaclib_db'
        character(len=*),parameter :: nuclib_db = 'nuclide_db'
        
        character(len=*), intent(in) :: datadir
        type(nuclib_data), intent(out) :: nuclib
        type(integer_dict), pointer :: nuclide_dict
        type(reaclib_data), intent(out) :: reaclib
        type(integer_dict), pointer :: rate_dict
        integer, intent(out) :: ierr

        character(len=160) :: reaclib_filename, nuclib_filename
        
        nuclib_filename = trim(datadir)//'/'//nuclib_db
        reaclib_filename = trim(datadir)//'/'//reaclib_db
        if (reaclib% Nentries /= 0) call free_reaclib_data(reaclib)

        ierr = 0
        write(error_unit,'(a)')  &
        & 'loading nuclib from '//trim(nuclib_filename)
        call do_load_nuclib(nuclib_filename,nuclib,ierr)
        write(error_unit,'(/,a,i0,a)')  &
        & 'done. ',nuclib% Nnuclides, &
        & ' nuclides retrieved. now writing nuclide dictionary...'
        call do_parse_nuclides(nuclib,nuclide_dict,ierr)
        write(error_unit,'(a//)') 'done.'
        
        write(error_unit,'(a)')  &
        & 'loading reaclib from '//trim(reaclib_filename)
        call do_load_reaclib(reaclib_filename,reaclib,ierr)
        write(error_unit,'(a,i0,a)') 'done. ',reaclib% Nentries, &
        & ' entries retrieved. now writing reaction dictionary...'
        call do_parse_rates(reaclib,rate_dict,ierr)
        write(error_unit,'(a,i0,a)') 'done. ',integer_dict_size(rate_dict), &
        & ' unique rates found.'
    end subroutine netJina_init

    subroutine get_nuclide_properties(nuclide, nuclib, nuclide_dict, &
            & A,Z,N,S,E,partition_fcn,provenance,ierr)
        use iso_fortran_env, only: error_unit
        use netJina_def
        use utils_lib, only: integer_dict_lookup

        character(len=iso_name_length), intent(in) :: nuclide
        type(nuclib_data), intent(in) :: nuclib
        type(integer_dict), pointer :: nuclide_dict
        real(dp), intent(out) :: A
        integer, intent(out) :: Z,N
        real(dp), intent(out) :: S,E
        real(dp), dimension(npfcn), intent(out) :: partition_fcn
        character(len=provenance_length), intent(out) :: provenance
        integer, intent(out) :: ierr
        integer :: indx
        
        A = 0.0
        Z = 0
        N = 0
        S = 0.0
        E = 0.0
        partition_fcn = 0.0
        provenance = ''
        
        call integer_dict_lookup(nuclide_dict,trim(nuclide),indx,ierr)
        if (ierr /= 0) then
            write(error_unit,'(a)') 'unable to find '//nuclide//' in database.'
            return
        end if
        
        ! debugging
        if (trim(nuclib% name(indx)) /= trim(nuclide)) then
            write(error_unit,'(a)') 'got the wrong nucleus'
            return
        end if
        
        provenance = nuclib% provenance(indx)
		A = nuclib% A(indx)
		Z = nuclib% Z(indx)
		N = nuclib% N(indx)
		S = nuclib% spin(indx)
		E = nuclib% mass_excess(indx)
		partition_fcn = nuclib% pfcn(:,indx)
    end subroutine get_nuclide_properties

    subroutine get_handle(reaclib,indx,handle)
        use netJina_def, only: reaclib_data, nJ_Nin, nJ_Nout, max_id_length
        use netJina_io, only: do_generate_handle
        
        type(reaclib_data), intent(in) :: reaclib
        integer, intent(in) :: indx
        character(len=max_id_length), intent(out) :: handle
        handle = do_generate_handle(reaclib% species(:,indx), &
        & nJ_Nin(reaclib% chapter(indx)), nJ_Nout(reaclib% chapter(indx)))
    end subroutine get_handle
    
    subroutine make_channel_handles(isotope,nuclib,nuclide_dict,handles,ierr)
        use netJina_def
        use utils_def
        use utils_lib
        use netJina_bdat
        character(len=iso_name_length), intent(in) :: isotope
        type(nuclib_data), intent(in) :: nuclib
        type(integer_dict), pointer :: nuclide_dict
        character(len=max_id_length),dimension(N_bdat_channels),intent(out) :: &
        & handles
        integer, intent(out) :: ierr
        call do_make_channel_handles(isotope,nuclib,nuclide_dict,handles,ierr)
    end subroutine make_channel_handles
    
    function reaction_string(reaclib,indx) result(str)
        use netJina_def
        type(reaclib_data), intent(in) :: reaclib
        integer, intent(in) :: indx
        character(len=length_reaction_string) :: str
        character(len=length_reaction_string) :: str_nxt
        
        str_nxt = ''
    	select case(reaclib% chapter(indx))
    		case(r_one_one)
    	   	call write_n_to_m(1,1)
    		case(r_one_two)
    	  	call write_n_to_m(1,2)
    		case(r_one_three)
    	  	call write_n_to_m(1,3)
    		case(r_two_one)
    	  	call write_n_to_m(2,1)
    		case(r_two_two)
    	  	call write_n_to_m(2,2)
    		case(r_two_three)
    	  	call write_n_to_m(2,3)
    		case(r_two_four)
    	  	call write_n_to_m(2,4)
    		case(r_three_one)
    	  	call write_n_to_m(3,1)
    		case(r_three_two)
    	  	call write_n_to_m(3,2)
    		case(r_four_two)
    	  	call write_n_to_m(4,2)
    		case(r_one_four)
    	  	call write_n_to_m(1,4)
    	end select
        
        str = str_nxt
    
    contains
    	subroutine write_n_to_m(n,m)
          	integer, intent(in) :: n, m
           	integer :: j
            do j = 1,n-1
                str_nxt = trim(str_nxt)//' '// &
                & trim(adjustl(reaclib% species(j, indx)))//' +'
        	end do
            str_nxt = trim(str_nxt)//' '// &
            & trim(adjustl(reaclib% species(n,indx)))//' ==>'
        	do j=n+1,n+m-1
        	    str_nxt = trim(str_nxt)//' '// &
        	    & trim(adjustl(reaclib% species(j, indx))) &
        	    & //' +' 
        	end do
            str_nxt = trim(str_nxt)//' '// &
            & trim(adjustl(reaclib% species(n+m,indx)))
    	end subroutine write_n_to_m	
        
    end function reaction_string
    
end module netJina_lib
