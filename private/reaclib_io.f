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
    
    subroutine do_read_reaclib(filename,rates,ierr)
    	use, intrinsic :: iso_fortran_env, only: iostat_end, error_unit
    	use netJina_def
        use utils_lib, only: alloc_iounit
        
    	character(len=*), intent(in) :: filename
    	type(reaclib_data), intent(out) :: rates
    	integer, intent(out) :: ierr
    	integer :: i, reaclib_unitno, count, iend

    	ierr = 0
        reaclib_unitno = alloc_iounit(ierr)
        if (ierr /= 0) return
        
    	open(unit=reaclib_unitno, file=trim(filename), iostat=ierr, status="old", action="read")
        if (failure('opening'//trim(filename),ierr)) return

    	! allocate a temporary to hold the library
    	call allocate_reaclib_data(rates,max_nreaclib,ierr)
        if (failure('allocating storage',ierr)) return

    	count = 0
    	do i = 1, max_nreaclib
    		read(unit=reaclib_unitno, fmt=reaclib_line0, iostat=iend) rates%chapter(i)
    		if (iend == iostat_end ) exit 
    		read(unit=reaclib_unitno,fmt=reaclib_line1,iostat=ierr) &
    		& rates%species(:,i), rates%label(i), rates%reaction_flag(i), &
    		& rates%reverse_flag(i),rates%Qvalue(i)
            if (failure('reading line 1/3',ierr)) return
    		
            read(unit=reaclib_unitno,fmt=reaclib_line2,iostat=ierr) &
    	    & rates%coefficients(1:4,i)
            if (failure('reading line 2/3',ierr)) return
    		
            read(unit=reaclib_unitno,fmt=reaclib_line3,iostat=ierr) &
    		& rates%coefficients(5:7,i)
            if (failure('reading line 3/3',ierr)) return
    		
    		count = count + 1
    	end do
    	close(reaclib_unitno)
        
    	rates% Nentries = count
        
    contains
        function failure(action,ierr)
            character(len=*), intent(in) :: action
            integer, intent(in) :: ierr
            logical :: failure
            character(len=*), parameter :: err_format = '("Error while ",a,i0)'
            
            if (ierr == 0) then
                failure = .TRUE.
                return
            end if
            failure = .FALSE.
            write (error_unit,err_format) action, ierr
        end function failure
    end subroutine do_read_reaclib
    

end module reaclib_io
