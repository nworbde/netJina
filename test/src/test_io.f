program test_io
    use iso_fortran_env, only: output_unit, error_unit
    use netJina_def
    use netJina_lib
    use utils_lib
    
    character(len=*), parameter :: datadir = '../data'
    type(reaclib_data) :: reaclib
    type(starlib_data) :: starlib
    type(nuclib_data) :: nuclib
    integer :: ierr
    integer :: i,j,nin,nout, Z, N
    real(dp) :: A, S, E, partition_fcn(npfcn)
    character(len=provenance_length) :: provenance
    type(integer_dict), pointer :: rates_dict=>null(), nuclide_dict=>null(), starlib_dict=>null()
    integer :: max_terms, indx(1), dindx, i_rate
    character(len=max_id_length) :: handle
    character(len=iso_name_length) :: iso
    character(len=max_id_length),dimension(N_bdat_channels) :: handles
    integer, dimension(N_bdat_channels) :: n_coeff
    real(dp), dimension(ncoefficients*max_terms_per_rate,N_bdat_channels) :: &
    & rate_coefficients
    real(dp), dimension(N_bdat_channels) :: q
    logical, dimension(N_bdat_channels) :: rate_mask
    real(dp), dimension(number_starlib_temps,N_bdat_channels) :: T9, rate, uncertainty
    
    call netJina_init(datadir,nuclib,nuclide_dict,reaclib, &
    &   starlib,rates_dict,starlib_dict,ierr)
    if (ierr /= 0) then
        write(error_unit,*) 'failure in initialization ',ierr
        stop
    end if

    write(output_unit,'(/,/,a)') 'What are the properties of cn337?'
    call get_nuclide_properties('cn337',nuclib,nuclide_dict, &
    & A,Z,N,S,E,partition_fcn,provenance,ierr)
    if (ierr /=0) then
        write(error_unit,*) 'cn337 is not on the list'
    else
        write(output_unit,'("ref.",t16,": ",a)') provenance
        write(output_unit,'("A, Z, N",t16,"=",f6.1,i4,i4)') A, Z, N
        write(output_unit,'("spin",t16,"=",f6.1)') S
        write(output_unit,'("mass excess",t16,"=",f10.3)') E
        write(output_unit,'(a,/,18("_"),/,a6,a12,/,18("="))')  &
        & 'partition function','T9','G(T9)'
        write(output_unit,'(f6.1,es12.4)') (pfcn_T9(i), partition_fcn(i),  &
        & i=1, npfcn)
        write(output_unit,'(18("_"))')
    end if
    
    do i = 1,reaclib% Nentries
        if (reaclib% chapter(i) == 5) exit
    end do
    if (reaclib% chapter(i) == 5) then
        write(output_unit,'(/,29("_"),/,a,/,29("="))')  &
        & 'First 50 entries of chapter 5'
        do j = i,i+50
            write(output_unit,'(a5,a5," ==> ",a5,a5)')  &
            & reaclib% species(1:4,j)
        end do
    else
        write(error_unit,'(a)') 'unable to find reactions in chapter 5'
    end if
    write(output_unit,'(29("_"))')

    write(output_unit,'(/,a)') 'What reaction has the most terms?'
    max_terms = maxval(reaclib% N_rate_terms)
    indx = maxloc(reaclib% N_rate_terms)
    write(output_unit,'(a,a,i3,a)') trim(reaction_string(reaclib,indx(1))), &
    & ' has ',max_terms,' terms:'
    call get_handle(reaclib,indx(1),handle)
    call integer_dict_lookup(rates_dict,trim(handle),dindx,ierr)
    if (ierr == 0) then
        do i = dindx-max_terms+1, dindx
            write(output_unit,'(i5,tr1,6a5,tr2,7es12.4)')  &
            & i,reaclib% species(:,i),reaclib% coefficients(:,i)
        end do
    end if

    write (output_unit,'(/,a)') 'What are the reaction channels for ca37?'
    iso = 'ca37'
    call make_channel_handles(iso,nuclib,nuclide_dict,handles,ierr)
    write(output_unit,'(54("_"),/,a2,tr2,a24,tr2,a24,/,54("="))')  &
    & 'id','handle','reaction'
    do i = 1, N_bdat_channels
        call integer_dict_lookup(rates_dict,trim(handles(i)),i_rate,ierr)
        if (ierr == 0) then
            write (output_unit,'(i2,tr2,a24,tr2,a24)') &
            &  i,trim(adjustl(handles(i))), &
            & trim(reaction_string(reaclib,i_rate))
        else
            write (output_unit,*) i,trim(handles(i)), 'rate not found'
        end if
    end do
    write(output_unit,'(54("_"))')

    write (output_unit,'(/,a)') 'with starlib...'
    iso = 'ca37'
    call make_channel_handles(iso,nuclib,nuclide_dict,handles,ierr)
    write(output_unit,'(54("_"),/,a2,tr2,a24,tr2,a24,/,54("="))')  &
    & 'id','handle','reaction'
    do i = 1, N_bdat_channels
        call integer_dict_lookup(starlib_dict,trim(handles(i)),i_rate,ierr)
        if (ierr == 0) then
            write (output_unit,'(i2,tr2,a24,tr2,a24)') &
            &  i,trim(adjustl(handles(i))), &
            & trim(reaction_string(starlib,i_rate))
        else
            write (output_unit,*) i,trim(handles(i)), 'rate not found'
        end if
    end do
    write(output_unit,'(54("_"))')

    write(output_unit,'(/,a)') 'What are the reaction parameters for ca37'
    call get_bdat_channels(reaclib,rates_dict, &
    & handles,n_coeff,rate_coefficients,q,rate_mask)
    write(output_unit,'(54("_"),/,a24,tr2,a10,tr2,a4,tr2,a7,/,54("="))') &
    & 'handle','q-value','mask','N terms'
    do i=1, N_bdat_channels
        write (output_unit,'(a24,tr2,f10.4,tr2,l4,tr2,i7)') &
        & trim(handles(i)),q(i),rate_mask(i),n_coeff(i)
    end do
    write(output_unit,'(54("_"))')

    write(output_unit,'(/,a)') 'with starlib...'
    call get_bdat_rates(starlib,starlib_dict, &
    & handles,T9,rate,uncertainty,q,rate_mask)
    write(output_unit,'(45("_"),/,a24,tr2,a10,tr2,a4,/,45("="))') &
    & 'handle','q-value','mask'
    do i=1, N_bdat_channels
        write (output_unit,'(a24,tr2,f10.4,tr2,l4)') &
        & trim(handles(i)),q(i),rate_mask(i)
    end do
    write(output_unit,'(54("_"))')
 
    write(output_unit,'(/,a)')  &
    & 'What are the reaction coefficients for ca37(g,p)k36?'
    write(output_unit,'(7(es12.4))') rate_coefficients(1:49,i_gp)

    write(output_unit,'(/,a)') &
    &   'What is the starlib rate for ca37(g,p)k36?'
    write(output_unit,'(2(es12.4))') (T9(j,i_gp),rate(j,i_gp), j=1,number_starlib_temps)
    
    write(output_unit,'(/,a)')  &
    & 'What are the reaction coefficients for ca37(n,a)ar34?'
    write(output_unit,'(7(es12.4))') rate_coefficients(1:7,i_na)

    write(output_unit,'(/,a)') &
    &   'What is the starlib rate for ca37(n,a)ar34?'
    write(output_unit,'(2(es12.4))') (T9(j,i_na),rate(j,i_na), j=1,number_starlib_temps)

end program test_io
