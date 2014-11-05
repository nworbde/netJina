!   reaclib_io.f
!
!   low-level read/write routines for the JINA reaclib database
!   Edward Brown, Michigan State University
!
!   Requires installation of MESA (mesa.sourceforge.net) utils lib, v 6480
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

    ! the starlib header format
    character(len=*), parameter :: header_format = '(a64)'
    
contains
    
    subroutine do_load_nuclib(filename,cache_filename,nuclides,ierr)
        use, intrinsic :: iso_fortran_env, only: iostat_end, error_unit
        use netJina_def
        use netJina_storage
        use utils_lib, only: alloc_iounit, free_iounit
        
        character(len=*), intent(in) :: filename,cache_filename
        type(nuclib_data), intent(out) :: nuclides
        integer, intent(out) :: ierr
        type(nuclib_data) :: tmp_n
        integer :: i, iso_count,ios, nuclib_unitno
        character(len=80) :: buffer
        real, dimension(24) :: pfcn
        logical :: have_cache
        
        ierr = 0
        inquire(file=cache_filename,exist=have_cache)
        if (have_cache) then
            call do_read_nuclib_cache(cache_filename,nuclides,ierr)
            if (ierr == 0) return
        end if
        
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
            tmp_n% name(i) = adjustl(tmp_n% name(i))
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
        
        if (.not.have_cache) then
            call do_write_nuclib_cache(cache_filename,nuclides,ierr)
        end if
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

    subroutine do_read_nuclib_cache(filename,nuclides,ierr)
        use netJina_def
        use netJina_storage
        character(len=*), intent(in) :: filename
        type(nuclib_data), intent(out) :: nuclides
        integer, intent(out) :: ierr
        integer :: cache_unitno, n
        
        ierr = 0
        
        open(newunit=cache_unitno, file=trim(filename), iostat=ierr, &
        &   action="read",status="old",form="unformatted")
        if (io_failure('opening '//trim(filename),ierr)) return
        
        read(cache_unitno) n
        call allocate_nuclib_data(nuclides,n, ierr)
        read(cache_unitno) nuclides% name
        read(cache_unitno) nuclides% provenance
        read(cache_unitno) nuclides% A
        read(cache_unitno) nuclides% Z
        read(cache_unitno) nuclides% N
        read(cache_unitno) nuclides% spin
        read(cache_unitno) nuclides% mass_excess
        read(cache_unitno) nuclides% pfcn
        read(cache_unitno) pfcn_T9
        close(cache_unitno)
    end subroutine do_read_nuclib_cache
     
    subroutine do_write_nuclib_cache(filename,nuclides,ierr)
        use netJina_def
        character(len=*), intent(in) :: filename
        type(nuclib_data), intent(in) :: nuclides
        integer, intent(out) :: ierr
        integer :: cache_unitno
        
        ierr = 0
        if (nuclides% Nnuclides == 0) return
        
        open(newunit=cache_unitno, file=trim(filename), iostat=ierr, &
        &   action="write",form="unformatted")
        if (io_failure('opening '//trim(filename),ierr)) return
        
        write(cache_unitno) nuclides% Nnuclides
        write(cache_unitno) nuclides% name
        write(cache_unitno) nuclides% provenance
        write(cache_unitno) nuclides% A
        write(cache_unitno) nuclides% Z
        write(cache_unitno) nuclides% N
        write(cache_unitno) nuclides% spin
        write(cache_unitno) nuclides% mass_excess
        write(cache_unitno) nuclides% pfcn
        write(cache_unitno) pfcn_T9

        close(cache_unitno)
    end subroutine do_write_nuclib_cache
   
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
    
    subroutine do_load_reaclib(filename,cache_filename,rates,ierr)
        use, intrinsic :: iso_fortran_env, only: iostat_end, error_unit
        use netJina_def
        use netJina_storage
        use utils_lib, only: alloc_iounit, free_iounit
        
        character(len=*), intent(in) :: filename, cache_filename
        type(reaclib_data), intent(out) :: rates
        integer, intent(out) :: ierr
        integer :: i, reaclib_unitno, count, iend
        type(reaclib_data) :: tmp_rates
        logical :: have_cache

        ierr = 0
        inquire(file=cache_filename,exist=have_cache)
        if (have_cache) then
            call do_read_reaclib_cache(cache_filename,rates,ierr)
            if (ierr == 0) return
        end if
        
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
            
            ! set flag if an electron capture rate
            if (trim(adjustl(tmp_rates% label(i))) == 'ec')  &
            &   tmp_rates% reaction_flag(i) = 'e'
            
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
        
        if (.not.have_cache) then
            call do_write_reaclib_cache(cache_filename,rates,ierr)
        end if
        
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
    
    subroutine do_read_reaclib_cache(filename,rates,ierr)
        use netJina_def
        use netJina_storage
        character(len=*), intent(in) :: filename
        type(reaclib_data), intent(out) :: rates
        integer, intent(out) :: ierr
        integer :: cache_unitno, n
        
        ierr = 0
        
        open(newunit=cache_unitno, file=trim(filename), iostat=ierr, &
        &   action="read",status="old",form="unformatted")
        if (io_failure('opening '//trim(filename),ierr)) return
        
        read(cache_unitno) n
        call allocate_reaclib_data(rates, n, ierr)
        read(cache_unitno) rates% chapter
        read(cache_unitno) rates% species
        read(cache_unitno) rates% label
        read(cache_unitno) rates% reaction_flag
        read(cache_unitno) rates% reverse_flag
        read(cache_unitno) rates% Qvalue
        read(cache_unitno) rates% coefficients
        read(cache_unitno) rates% N_rate_terms

        close(cache_unitno)
    end subroutine do_read_reaclib_cache
     
    subroutine do_write_reaclib_cache(filename,rates,ierr)
        use netJina_def
        character(len=*), intent(in) :: filename
        type(reaclib_data), intent(in) :: rates
        integer, intent(out) :: ierr
        integer :: cache_unitno
        
        ierr = 0
        if (rates% Nentries == 0) return
        
        open(newunit=cache_unitno, file=trim(filename), iostat=ierr, &
        &   action="write",form="unformatted")
        if (io_failure('opening '//trim(filename),ierr)) return
        
        write(cache_unitno) rates% Nentries
        write(cache_unitno) rates% chapter
        write(cache_unitno) rates% species
        write(cache_unitno) rates% label
        write(cache_unitno) rates% reaction_flag
        write(cache_unitno) rates% reverse_flag
        write(cache_unitno) rates% Qvalue
        write(cache_unitno) rates% coefficients
        write(cache_unitno) rates% N_rate_terms

        close(cache_unitno)
    end subroutine do_write_reaclib_cache
    
    subroutine do_parse_rates(reaclib,rate_dict,ierr)
        use, intrinsic :: iso_fortran_env, only : error_unit
        use utils_def, only: integer_dict
        use utils_lib, only: integer_dict_define_and_report_duplicates, integer_dict_free
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
        ierr = 0
        do i_rate = 1, reaclib% Nentries
            nin = nJ_Nin(reaclib% chapter(i_rate))
            nout = nJ_Nout(reaclib% chapter(i_rate))
            handle = do_generate_handle(reaclib% species(:,i_rate), nin, nout)

            ! p+p->d is special
            if (trim(handle) == 'nJ_p_p_to_d' .and.  &
            &   reaclib% reaction_flag(i_rate) == 'e') then
                handle = trim(handle)//'_e'
            end if

            call integer_dict_define_and_report_duplicates(rate_dict,trim(handle),i_rate,duplicate,ierr)
            if (.not. duplicate) then
                reaclib% N_rate_terms(i_rate) = 1
            else
                reaclib% N_rate_terms(i_rate) = reaclib% N_rate_terms(i_rate-1) + 1
            end if
            if (mod(i_rate,1000) == 0)  &
            & write (error_unit,'(a)',advance='no') '.'
        end do
    end subroutine do_parse_rates
    
    subroutine do_load_starlib(filename,cache_filename,rates,ierr)
        use, intrinsic :: iso_fortran_env, only: iostat_end, error_unit
        use netJina_def
        use netJina_storage
        use utils_lib, only: alloc_iounit, free_iounit
        
        character(len=*), intent(in) :: filename,cache_filename
        type(starlib_data), intent(out) :: rates
        integer, intent(out) :: ierr
        integer :: i, starlib_unitno, count, iend, j
        type(starlib_data) :: tmp_rates
        character(len=64) :: header
        real(dp), dimension(number_starlib_temps) :: T9,rate,uncertainty
        integer :: chapter
        integer :: ns
        character(len=iso_name_length),  &
        &   dimension(starlib_max_species_per_reaction) :: &
        & species
        character(len=reaction_reference_length) :: prov
        character :: reverse, rflag
        real(dp) :: q
        logical :: have_cache
        
        ierr = 0
        inquire(file=cache_filename,exist=have_cache)
        if (have_cache) then
            call do_read_starlib_cache(cache_filename,rates,ierr)
            if (ierr == 0) return
        end if

        starlib_unitno = alloc_iounit(ierr)
        if (io_failure('allocating iounit',ierr)) return
        
        open(unit=starlib_unitno, file=trim(filename), iostat=ierr, &
        & status="old", action="read")
        if (io_failure('opening '//trim(filename),ierr)) return

        ! allocate a temporary to hold the library
        call allocate_starlib_data(tmp_rates,max_nstarlib,ierr)
        if (io_failure('allocating storage',ierr)) return

        count = 0
        do i = 1, max_nstarlib
            read(unit=starlib_unitno, fmt=header_format, iostat=iend) header
            if (iend == iostat_end ) exit
            if (io_failure('reading starlib header',ierr)) return
            call parse_header(header,chapter,ns,species,prov,rflag,reverse,q)
            
            read(unit=starlib_unitno,fmt=*,iostat=ierr) &
            &    (T9(j),rate(j),uncertainty(j),j=1,number_starlib_temps)
            if (io_failure('reading starlib rates',ierr)) return
            
            tmp_rates% chapter(i) = chapter
            tmp_rates% species(:,i) = species
            tmp_rates% label(i) = prov
            tmp_rates% reaction_flag(i) = rflag
            tmp_rates% reverse_flag(i) = reverse
            tmp_rates% Qvalue(i) = q
            tmp_rates% T9(:,i) = T9
            tmp_rates% rate(:,i) = rate
            tmp_rates% uncertainty(:,i) = uncertainty
            count = count + 1
            if (mod(count,1000) == 0)  &
            &   write (error_unit,'(a)',advance='no') '.'
            
        end do
        close(starlib_unitno)
        call free_iounit(starlib_unitno)
        tmp_rates% Nentries = count
                
        call copy_starlib_data(tmp_rates,rates,ierr)
        if (io_failure('copying starlib data',ierr)) return
        call free_starlib_data(tmp_rates)
        
        if (.not.have_cache) then
            call do_write_starlib_cache(cache_filename,rates,ierr)
        end if
        
    contains
        subroutine parse_header(header,chapter,ns,species,prov,rflag,reverse,q)
            character(len=64), intent(in) :: header
            integer, intent(out) :: chapter
            integer, intent(out) :: ns
            character(len=iso_name_length), &
            &  dimension(starlib_max_species_per_reaction), intent(out) :: &
            &    species
            character(len=reaction_reference_length), intent(out) :: prov
            character, intent(out) :: rflag, reverse
            real(dp), intent(out) :: q
            character(len=iso_name_length), &
            &   dimension(starlib_max_species_per_reaction) :: tmp_species
            integer :: j,i
        
            read(header,'(i5,2(a5),4(a5),t44,a4,a1,t53,es12.5)')  &
            &   chapter, tmp_species(:), prov, rflag, q
        
            ! check if reaction flag indicates a reverse rate
            reverse = ''
            if (rflag == 'v') reverse = rflag
            
            ! check for electron captures
            if (trim(adjustl(prov)) == 'ec') rflag = 'e'
            
            species = ''
            ns = 0
            do j= 1, starlib_max_species_per_reaction
                if (len_trim(tmp_species(j)) > 0) then
                    ns = ns + 1
                    species(ns) = adjustl(tmp_species(j))
                end if
            end do
        
        end subroutine parse_header
     
        subroutine copy_starlib_data(old_data,new_data,ierr)
            type(starlib_data), intent(in) :: old_data
            type(starlib_data), intent(out) :: new_data
            integer, intent(out) :: ierr
            integer :: n
        
            if (new_data% Nentries /= 0) call free_starlib_data(new_data)
            call allocate_starlib_data(new_data,old_data% Nentries, ierr)
            if (ierr /= 0) return
        
            n = old_data% Nentries
            new_data% chapter(:) = old_data% chapter(1:n)
            new_data% species(:,:) = old_data% species(:,1:n)
            new_data% label(:) = old_data% label(1:n)
            new_data% reaction_flag(:) = old_data% reaction_flag(1:n)
            new_data% reverse_flag(:) = old_data% reverse_flag(1:n)
            new_data% Qvalue(:) = old_data% Qvalue(1:n)
            new_data% T9(:,:) = old_data% T9(:,1:n)
            new_data% rate(:,:) = old_data% rate(:,1:n)
            new_data% uncertainty(:,:) = old_data% uncertainty(:,1:n)
        end subroutine copy_starlib_data
    end subroutine do_load_starlib    
    
    subroutine do_read_starlib_cache(filename,rates,ierr)
        use netJina_def
        use netJina_storage
        character(len=*), intent(in) :: filename
        type(starlib_data), intent(out) :: rates
        integer, intent(out) :: ierr
        integer :: cache_unitno, n
        
        ierr = 0
        
        open(newunit=cache_unitno, file=trim(filename), iostat=ierr, &
        &   action="read",status="old",form="unformatted")
        if (io_failure('opening '//trim(filename),ierr)) return
        
        read(cache_unitno) n
        call allocate_starlib_data(rates, n, ierr)
        read(cache_unitno) rates% chapter
        read(cache_unitno) rates% species
        read(cache_unitno) rates% label
        read(cache_unitno) rates% reaction_flag
        read(cache_unitno) rates% reverse_flag
        read(cache_unitno) rates% Qvalue
        read(cache_unitno) rates% T9
        read(cache_unitno) rates% rate
        read(cache_unitno) rates% uncertainty

        close(cache_unitno)
    end subroutine do_read_starlib_cache
     
    subroutine do_write_starlib_cache(filename,rates,ierr)
        use netJina_def
        character(len=*), intent(in) :: filename
        type(starlib_data), intent(in) :: rates
        integer, intent(out) :: ierr
        integer :: cache_unitno
        
        ierr = 0
        if (rates% Nentries == 0) return
        
        open(newunit=cache_unitno, file=trim(filename), iostat=ierr, &
        &   action="write",form="unformatted")
        if (io_failure('opening '//trim(filename),ierr)) return
        
        write(cache_unitno) rates% Nentries
        write(cache_unitno) rates% chapter
        write(cache_unitno) rates% species
        write(cache_unitno) rates% label
        write(cache_unitno) rates% reaction_flag
        write(cache_unitno) rates% reverse_flag
        write(cache_unitno) rates% Qvalue
        write(cache_unitno) rates% T9
        write(cache_unitno) rates% rate
        write(cache_unitno) rates% uncertainty

        close(cache_unitno)
    end subroutine do_write_starlib_cache
    
    subroutine do_parse_starlib(starlib,starlib_dict,ierr)
        use, intrinsic :: iso_fortran_env, only : error_unit
        use utils_def, only: integer_dict
        use utils_lib, only: integer_dict_define_and_report_duplicates, integer_dict_free
        use netJina_def
        
        type(starlib_data), intent(inout) :: starlib
        type(integer_dict), pointer:: starlib_dict
        integer, intent(out) :: ierr
        integer :: i_rate, nin, nout, ikey, indx, nterms
        character(len=max_id_length) :: handle
        integer :: i_head
        logical :: duplicate
        
        if (associated(starlib_dict)) call integer_dict_free(starlib_dict)      
        ierr = 0
        do i_rate = 1, starlib% Nentries
            nin = nJ_Nin(starlib% chapter(i_rate))
            nout = nJ_Nout(starlib% chapter(i_rate))
            handle = do_generate_handle(starlib% species(:,i_rate), nin, nout)
 
            ! p+p->d is special
            if (trim(handle) == 'nJ_p_p_to_d' .and.  &
            &   starlib% reaction_flag(i_rate) == 'e') then
                handle = trim(handle)//'_e'
            end if
 
            call integer_dict_define_and_report_duplicates(starlib_dict, &
            &   trim(handle),i_rate,duplicate,ierr)
            if (duplicate) then
                write(error_unit,'(a)') trim(handle)//' is duplicate'
            end if
            if (mod(i_rate,1000) == 0)  &
            & write (error_unit,'(a)',advance='no') '.'
        end do
    end subroutine do_parse_starlib
    
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
