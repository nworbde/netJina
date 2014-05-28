!   reaclib_io.f
!
!   low-level read/write routines for the JINA reaclib database
!   Edward Brown, Michigan State University
!
!   Requires installation of MESA (mesa.sourceforge.net) utils lib, v 6022
!

module netJina_io

    ! the winvne format
    character(len=*), parameter :: pfcn_fmt = '(23f3.2,f3.1)'
    character(len=*), parameter :: iso_fmt = '(a5)'
    character(len=*), parameter :: line0_fmt = '(a5,f12.3,2i4,f6.1,f10.3,a6)'

    ! the reaclib v2 data format
    character(len=*), parameter :: reaclib_line0 = '(i2)'
    character(len=*), parameter :: reaclib_line1 = &
    & '(5x,6a5,8x,a4,a1,a1,3x,es12.5)'
    character(len=*), parameter :: reaclib_line2 = '(4es13.6)'
    character(len=*), parameter :: reaclib_line3 = '(3es13.6)'

contains
    
    subroutine do_load_nuclib(filename,nuclides,ierr)
        use, intrinsic :: iso_fortran_env, only: iostat_end, error_unit
        use netJina_def
        use netJina_storage
        use utils_lib, only: alloc_iounit, free_iounit
        
        character(len=*), intent(in) :: filename
        type(nuclib_data), intent(out) :: nuclides
        integer, intent(out) :: ierr
        type(nuclib_data) :: tmp_n
        integer :: i, iso_count,ios, nuclib_unitno
        character(len=80) :: buffer
        real, dimension(24) :: pfcn
    
        ierr = 0
        nuclib_unitno = alloc_iounit(ierr)
        if (io_failure('allocating iounit',ierr)) return
        
        open(unit=nuclib_unitno, file=trim(filename), iostat=ierr, &
        & status="old", action="read")
        if (io_failure('opening '//trim(filename),ierr)) return

        ! allocate a temporary to hold the library
        call allocate_nuclib_data(tmp_n,max_nnuclib,ierr)
        if (io_failure('allocating storage',ierr)) return
    
        ! read in the partition function temperatures
        read(nuclib_unitno,pfcn_fmt,iostat=ierr) pfcn_T9(:)
        if (io_failure('reading partition function',ierr)) return

        i = 1
        do
            ! read in the full data line and check for length
            read(nuclib_unitno,fmt='(a80)',iostat=ios) buffer
            if (ios == iostat_end) exit
            if (len_trim(buffer) <= iso_name_length) cycle
            read(buffer,line0_fmt) &
            & tmp_n% name(i), tmp_n% A(i), tmp_n% Z(i), tmp_n% N(i),  &
            & tmp_n% spin(i), tmp_n% mass_excess(i),tmp_n% provenance(i)
            read(nuclib_unitno,*) tmp_n% pfcn(1:8,i)
            read(nuclib_unitno,*) tmp_n% pfcn(9:16,i)
            read(nuclib_unitno,*) tmp_n% pfcn(17:24,i)
            i = i+1
            if (mod(i,100) == 0)  &
            & write (error_unit,'(a)',advance='no') '.'
        end do
        close(nuclib_unitno)
        call free_iounit(nuclib_unitno)
        tmp_n% Nnuclides = i-1
        call copy_nuclib_data(tmp_n,nuclides,ierr)
        if (io_failure('copying nuclib data',ierr)) return
        call free_nuclib_data(tmp_n)
    contains
        
        subroutine copy_nuclib_data(old_data,new_data,ierr)
            type(nuclib_data), intent(in) :: old_data
            type(nuclib_data), intent(out) :: new_data
            integer, intent(out) :: ierr
            integer :: n
        
            if (new_data% Nnuclides /= 0) call free_nuclib_data(new_data)
            call allocate_nuclib_data(new_data,old_data% Nnuclides, ierr)
            if (ierr /= 0) return
        
            n = old_data% Nnuclides
            new_data% name(:) = old_data% name(1:n)
            new_data% provenance(:) = old_data% provenance(1:n)
            new_data% A(:) = old_data% A(1:n)
            new_data% Z(:) = old_data% Z(1:n)
            new_data% N(:) = old_data% N(1:n)
            new_data% spin(:) = old_data% spin(1:n)
            new_data% mass_excess(:) = old_data% mass_excess(1:n)
            new_data% pfcn(:,:) = old_data% pfcn(:,1:n)
        end subroutine copy_nuclib_data
        
    end subroutine do_load_nuclib    
    
    subroutine do_parse_nuclides(nuclib,nuclide_dict,ierr)
        use utils_def, only: integer_dict
        use utils_lib, only: integer_dict_define, &
        & integer_dict_lookup, integer_dict_free, integer_dict_create_hash
        use netJina_def
        
        type(nuclib_data), intent(in) :: nuclib
        type(integer_dict), pointer:: nuclide_dict
        integer, intent(out) :: ierr
        integer :: indx, ikey
         
        if (associated(nuclide_dict)) call integer_dict_free(nuclide_dict)
        do indx = 1, nuclib% Nnuclides
            call integer_dict_define(nuclide_dict, &
            & trim(adjustl(nuclib% name(indx))),indx,ierr)
        end do
        call integer_dict_create_hash(nuclide_dict,ierr)
    end subroutine do_parse_nuclides
    
    subroutine do_load_reaclib(filename,rates,ierr)
        use, intrinsic :: iso_fortran_env, only: iostat_end, error_unit
        use netJina_def
        use netJina_storage
        use utils_lib, only: alloc_iounit, free_iounit
        
        character(len=*), intent(in) :: filename
        type(reaclib_data), intent(out) :: rates
        integer, intent(out) :: ierr
        integer :: i, reaclib_unitno, count, iend
        type(reaclib_data) :: tmp_rates

        ierr = 0
        reaclib_unitno = alloc_iounit(ierr)
        if (io_failure('allocating iounit',ierr)) return
        
        open(unit=reaclib_unitno, file=trim(filename), iostat=ierr, &
        & status="old", action="read")
        if (io_failure('opening '//trim(filename),ierr)) return

        ! allocate a temporary to hold the library
        call allocate_reaclib_data(tmp_rates,max_nreaclib,ierr)
        if (io_failure('allocating storage',ierr)) return

        count = 0
        do i = 1, max_nreaclib
            read(unit=reaclib_unitno, fmt=reaclib_line0, iostat=iend) &
            & tmp_rates% chapter(i)
            if (iend == iostat_end ) exit 
            read(unit=reaclib_unitno,fmt=reaclib_line1,iostat=ierr) &
            & tmp_rates% species(:,i), tmp_rates% label(i),  &
            & tmp_rates% reaction_flag(i), &
            & tmp_rates% reverse_flag(i),tmp_rates% Qvalue(i)
            if (io_failure('reading line 1/3',ierr)) return
            
            read(unit=reaclib_unitno,fmt=reaclib_line2,iostat=ierr) &
            & tmp_rates% coefficients(1:4,i)
            if (io_failure('reading line 2/3',ierr)) return
            
            read(unit=reaclib_unitno,fmt=reaclib_line3,iostat=ierr) &
            & tmp_rates% coefficients(5:7,i)
            if (io_failure('reading line 3/3',ierr)) return
            
            count = count + 1
            if (mod(count,1000) == 0)  &
            & write (error_unit,'(a)',advance='no') '.'
            
        end do
        close(reaclib_unitno)
        call free_iounit(reaclib_unitno)
        tmp_rates%  Nentries = count
                
        call copy_reaclib_data(tmp_rates,rates,ierr)
        if (io_failure('copying reaclib data',ierr)) return
        
        call free_reaclib_data(tmp_rates)
        
    contains
     
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
        use, intrinsic :: iso_fortran_env, only : error_unit
        use utils_def, only: integer_dict
        use utils_lib, only: integer_dict_define_and_check
        use netJina_def
        
        type(reaclib_data), intent(inout) :: reaclib
        type(integer_dict), pointer:: rate_dict
        integer, intent(out) :: ierr
        integer :: i_rate, nin, nout, ikey, indx, nterms
        character(len=max_id_length) :: handle
        integer :: i_head
        logical :: duplicate
        
        if (associated(rate_dict)) call integer_dict_free(rate_dict)
        reaclib% N_rate_terms(:) = 1
        write (error_unit, '(a)') 'generating handles'        

        do i_rate = 1, reaclib% Nentries
            nin = nJ_Nin(reaclib% chapter(i_rate))
            nout = nJ_Nout(reaclib% chapter(i_rate))
            handle = do_generate_handle(reaclib% species(:,i_rate), nin, nout)

            call integer_dict_define_and_check(rate_dict,trim(handle),i_rate,duplicate,ierr)
            if (.not. duplicate) then
                reaclib% N_rate_terms(i_rate) = 1
            else
                reaclib% N_rate_terms(i_rate) = reaclib% N_rate_terms(i_rate-1) + 1
            end if
            if (mod(i_rate,1000) == 0)  &
            & write (error_unit,'(a)',advance='no') '.'
        end do
        write (error_unit,'(a)') 'done'
    end subroutine do_parse_rates
    
    function do_generate_handle(species,nin,nout) result(handle)
        use netJina_def
        
        character(len=iso_name_length),dimension(:),intent(in) :: species
        integer, intent(in) :: nin, nout
        character(len=max_id_length) :: handle
        character(len=2),parameter :: nj = 'nJ'
        character(len=3),parameter :: sep = '_to'
        integer :: i
        handle = nj
        do i = 1,nin
            handle = trim(handle)//'_'//adjustl(species(i))
        end do
        handle = trim(handle)//sep
        do i = nin+1,min(nin+nout,size(species))
            handle = trim(handle)//'_'//adjustl(species(i))
        end do
    end function do_generate_handle
    
    function io_failure(action,ierr)
        use iso_fortran_env, only: error_unit
        character(len=*), intent(in) :: action
        integer, intent(in) :: ierr
        logical :: io_failure
        character(len=*), parameter :: err_format = &
        &  '("Error while ",a,". error code = ",i0)'
        
        if (ierr == 0) then
            io_failure = .FALSE.
            return
        end if
        io_failure = .TRUE.
        write (error_unit,err_format) action, ierr
    end function io_failure
    
end module netJina_io
