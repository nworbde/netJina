!   reaclib_io.f
!
!   low-level read/write routines for the JINA reaclib database
!   Edward Brown, Michigan State University
!
!   Requires installation of MESA (mesa.sourceforge.net) utils lib, v 6022
!

module reaclib_io

    ! the reaclib v2 data format
	character(len=*), parameter :: reaclib_line0 = '(i2)'
	character(len=*), parameter :: reaclib_line1 = &
	& '(5x,6a5,8x,a4,a1,a1,3x,es12.5)'
	character(len=*), parameter :: reaclib_line2 = '(4es13.6)'
	character(len=*), parameter :: reaclib_line3 = '(3es13.6)'

contains
    
    subroutine do_load_reaclib(filename,rates,ierr)
    	use, intrinsic :: iso_fortran_env, only: iostat_end, error_unit
    	use netJina_def
        use utils_lib, only: alloc_iounit
        
    	character(len=*), intent(in) :: filename
    	type(reaclib_data), intent(out) :: rates
    	integer, intent(out) :: ierr
    	integer :: i, reaclib_unitno, count, iend
        type(reaclib_data) :: tmp_rates

    	ierr = 0
        reaclib_unitno = alloc_iounit(ierr)
        if (failure('allocating iounit',ierr)) return
        
    	open(unit=reaclib_unitno, file=trim(filename), iostat=ierr, &
    	& status="old", action="read")
        if (failure('opening'//trim(filename),ierr)) return

    	! allocate a temporary to hold the library
    	call allocate_reaclib_data(tmp_rates,max_nreaclib,ierr)
        if (failure('allocating storage',ierr)) return

    	count = 0
    	do i = 1, max_nreaclib
    		read(unit=reaclib_unitno, fmt=reaclib_line0, iostat=iend) &
    		& tmp_rates% chapter(i)
    		if (iend == iostat_end ) exit 
    		read(unit=reaclib_unitno,fmt=reaclib_line1,iostat=ierr) &
    		& tmp_rates% species(:,i), tmp_rates% label(i),  &
    		& tmp_rates% reaction_flag(i), &
    		& tmp_rates% reverse_flag(i),tmp_rates% Qvalue(i)
            if (failure('reading line 1/3',ierr)) return
    		
            read(unit=reaclib_unitno,fmt=reaclib_line2,iostat=ierr) &
    	    & tmp_rates% coefficients(1:4,i)
            if (failure('reading line 2/3',ierr)) return
    		
            read(unit=reaclib_unitno,fmt=reaclib_line3,iostat=ierr) &
    		& tmp_rates% coefficients(5:7,i)
            if (failure('reading line 3/3',ierr)) return
    		
    		count = count + 1
    	end do
    	close(reaclib_unitno)
        write(error_unit,*) 'received ',count,' entries'
    	tmp_rates%  Nentries = count
                
        call copy_reaclib_data(tmp_rates,rates,ierr)
        if (failure('copying reaclib data',ierr)) return
        
        call free_reaclib_data(tmp_rates)
        
    contains
        function failure(action,ierr)
            character(len=*), intent(in) :: action
            integer, intent(in) :: ierr
            logical :: failure
            character(len=*), parameter :: err_format = '("Error while ",a,i0)'
            
            if (ierr == 0) then
                failure = .FALSE.
                return
            end if
            failure = .TRUE.
            write (error_unit,err_format) action, ierr
        end function failure
        
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
    end subroutine do_load_reaclib
    
    subroutine do_parse_rates(reaclib,rate_dict,ierr)
        use utils_def, only: integer_dict
        use utils_lib, only: integer_dict_define, &
        & integer_dict_lookup, integer_dict_free
        use netJina_def
        
        type(reaclib_data), intent(inout) :: reaclib
        type(integer_dict), pointer:: rate_dict
        integer, intent(out) :: ierr
        integer :: i_rate, nin, nout, ikey, indx, nterms
        character(len=max_id_length) :: handle
         
        if (associated(rate_dict)) call integer_dict_free(rate_dict)
        reaclib% N_rate_terms(:) = 1
        do i_rate = 1, reaclib% Nentries
            nin = nJ_Nin(reaclib% chapter(i_rate))
            nout = nJ_Nout(reaclib% chapter(i_rate))
            handle = generate_handle(reaclib% species(:,i_rate), Nin, Nout)
            ! if we have a new handle, enter the location in the dictionary
            ! otherwise, increment the N_rate_terms entry
            call integer_dict_lookup(rate_dict,handle,indx,ikey)
            if (ikey /= 0) then
                call integer_dict_define(rate_dict,handle,i_rate,ierr)
            else
                reaclib% N_rate_terms(indx) = reaclib% N_rate_terms(indx) + 1
            end if
        end do
        
        ! now pass through and set N_rate_terms for all terms
        i_rate = 1
        do
            if (i_rate >= reaclib% Nentries) exit
            nterms = reaclib% N_rate_terms(i_rate)
            if (reaclib% N_rate_terms(i_rate) > 1) then
                reaclib% N_rate_terms(i_rate+1:i_rate+nterms-1) =  &
                & reaclib% N_rate_terms(i_rate)
            end if
            i_rate = i_rate+nterms
        end do
    end subroutine do_parse_rates
    
    function generate_handle(species,nin,nout) result(handle)
        use netJina_def
        
        character(len=iso_name_length),dimension(:),intent(in) :: species
        integer, intent(in) :: nin, nout
        character(len=max_id_length) :: handle
        character(len=2),parameter :: nj = 'nJ'
        character(len=4),parameter :: sep = '_to_'
        integer :: i
        handle = nj
        do i = 1,nin
            handle = trim(handle)//adjustl(species(i))
        end do
        handle = trim(handle)//sep
        do i = nin+1,min(nin+nout,size(species))
            handle = trim(handle)//adjustl(species(i))
        end do
    end function generate_handle
    
end module reaclib_io
