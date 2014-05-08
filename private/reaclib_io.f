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
    end subroutine do_load_reaclib
    

end module reaclib_io
