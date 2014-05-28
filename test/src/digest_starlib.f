program digest_starlib
    use, intrinsic :: iso_fortran_env, only: error_unit, iostat_end
    integer, parameter :: max_number_rates = 80000
    integer, parameter :: number_temps = 60
    character(len=*), parameter :: header_format = '(a64)'
    character(len=64) :: header
    real(8), dimension(number_temps) :: T9,rate,uncertainty
    real(8) :: foo
    integer :: nrates,iounit,ios, i, j
    integer :: ns
    character(len=5),dimension(4) :: species
    character(len=4) :: prov
    logical :: reverse
    real(8) :: q
    
    open(newunit=iounit,file='../data/starlib_db', &
    &   status='old',action='read',iostat=ios)
    if (ios /= 0) then
        write(error_unit,*) 'unable to open starlib_db'
        stop
    end if
    
    nrates = 0
    do i = 1, max_number_rates
        read (iounit,header_format,iostat=ios) header
        if (ios == iostat_end) then
            exit
        end if
        if (ios /= 0) then
            write(error_unit,*) 'error in header read at rate ',nrates
            stop
        end if
        read(iounit,*) (T9(j),rate(j),uncertainty(j),j=1,number_temps)
        if (mod(nrates,1000) == 0) then
            write(error_unit,*) header
            call parse_header(header,ns,species,prov,reverse,q)
            write (error_unit,'(i6,tr2,i2,tr2,4(a5,tr1),a4,tr1,l1,tr1,es11.4)') i,ns,species,prov,reverse,q
            write(error_unit,*) T9(1),rate(1)
        end if
        nrates = nrates + 1
    end do

    close(iounit)
    write (*,*) 'number of rates = ',nrates
    
contains
    subroutine parse_header(header,ns,species,prov,reverse,q)
        character(len=64), intent(in) :: header
        integer, intent(out) :: ns
        character(len=5),dimension(4), intent(out) :: species
        character(len=4), intent(out) :: prov
        logical, intent(out) :: reverse
        real(8), intent(out) :: q
        integer :: chapter
        character :: rf
        character(len=5),dimension(6) :: tmp_species
        integer :: j,i
        
        read(header,'(i5,2(a5),4(a5),t44,a4,a1,t53,es12.5)')  &
        &   chapter, tmp_species(:), prov, rf, q
        
        reverse = .FALSE.
        if (rf == 'v') reverse = .TRUE.

        species = ''
        ns = 0
        do j= 1, 6
            if (len_trim(tmp_species(j)) > 0) then
                ns = ns + 1
                species(ns) = adjustl(tmp_species(j))
            end if
        end do
        
    end subroutine parse_header
end program digest_starlib
