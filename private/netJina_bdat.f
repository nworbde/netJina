module netJina_bdat

contains
    subroutine make_channel_handles(isotope,nuclib,nuclide_dict,handles,ierr)
        use netJina_def
        use utils_def
        use utils_lib
        character(len=iso_name_length), intent(in) :: isotope,
        type(nuclib_data), intent(in) :: nuclib,
        type(integer_dict), pointer :: nuclide_dict
        charactet(len=max_id_length),dimension(N_bdat_channels),intent(out) :: &
        & handles
        integer, intent(out) :: ierr
        integer :: indx
        character(len=iso_name_length), dimension(N_bdat_channels) :: product
        integer :: Z, N, Zt, Nt
        
        ! lookup nuclide
        callinteger_dict_lookup(nuclide_dict,isotope,indx,ierr)
        if (failure('unable to find '//trim(isotope)//'. error code ',ierr)) &
        & return
        Z = nuclib% Z(indx)
        N = nuclib% N(indx)
        
        ! (p,g)
        Zt = Z + 1
        Nt = N
        write(product(i_pg),'(a0,i0)') element_name(Zt),Zt+Nt
        ! (a,n)
        Zt = Z + 2
        Nt = N + 1
        write(product(i_an),'(a0,i0)') element_name(Zt),Zt+Nt
        ! (a,g)
        Zt = Z + 2
        Nt = N + 2
        write(product(i_ag),'(a0,i0)') element_name(Zt),Zt+Nt
        ! (a,p)
        Zt = Z + 1
        Nt = N + 2
        write(product(i_ap),'(a0,i0)') element_name(Zt),Zt+Nt
        ! (n,g)
        Zt = Z
        Nt = N + 1
        write(product(i_ng),'(a0,i0)') element_name(Zt),Zt+Nt
        ! (n,p)
        Zt = Z - 1
        Nt = N + 1
        write(product(i_np),'(a0,i0)') element_name(Zt),Zt+Nt
        ! (g,p)
        Zt = Z - 1
        Nt = N
        write(product(i_gp),'(a0,i0)') element_name(Zt),Zt+Nt
        ! (n,a)
        Zt = Z - 2
        Nt = N - 1
        write(product(i_na),'(a0,i0)') element_name(Zt),Zt+Nt
        ! (g,a)
        Zt = Z - 2
        Nt = N - 2
        write(product(i_ga),'(a0,i0)') element_name(Zt),Zt+Nt
        ! (p,a)
        Zt = Z - 1
        Nt = N - 2
        write(product(i_pa),'(a0,i0)') element_name(Zt),Zt+Nt
        ! (g,n)
        Zt = Z
        Nt = N - 1
        write(product(i_gn),'(a0,i0)') element_name(Zt),Zt+Nt
        ! (p,n)
        Zt = Z + 1
        Nt = N - 1
        write(product(i_pn),'(a0,i0)') element_name(Zt),Zt+Nt
        ! (p,g)
        Zt = Z + 1
        Nt = N
        write(product(i_pg),'(a0,i0)') element_name(Zt),Zt+Nt
        
    end subroutine make_channel_handles

    

end module netJina_bdat
