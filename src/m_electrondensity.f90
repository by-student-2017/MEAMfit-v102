module m_electrondensity !f_i_L_ij,rho_i_l
    !f_j_L_ij, computed in subroutine radialdensityfunction
    real(8), allocatable:: fjlij(:,:,:),   & !L,site j,site i
        rhol(:,:),      & !site,L
        rho_i(:)          !site
end module m_electrondensity
