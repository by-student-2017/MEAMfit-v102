module m_meamparameters
    !MEAM parameters
    logical cmin_cmax_zero,lookuptables,pairpotonly,envdepmeamt
    integer lmax!,nCutoffOverride
    real(8), allocatable:: cmin(:,:,:), & !species 3x  [,]
        cmax(:,:,:),         & !species 3x  [,]
        meamtau(:,:),        & !species,L [,]
        meamrhodecay(:,:,:,:),& !6,L,species,species
        meame0(:),           & !species  [,]
        meamrho0(:),         & !species  [,]
        meamemb3(:),           & !species  [,]
        meamemb4(:),         & !species  [,]
        rs(:,:),             & !species 2x
        rc(:,:),             & !species 2x
        !temporary values &
    meam_f(:),           & !site
        meam_t(:,:),         & !site,L
        meam_paire(:),       & !site
        enconst(:)           !species
    !Pair potential parameters
    integer, parameter:: pot_samespecies_arrsize=1000 !For internal
    !Fe-Fe Ackland potential 1000
    integer emb_Fe_arrsz,emb_C_arrsz !For external look-up tables
    real(8), allocatable:: pairpotparameter(:,:,:) !ispecies,jspecies
    real(8), allocatable:: pot_samespecies(:),sec_der(:) !For internal
    ! Fe-Fe Ackland potential
    real(8) r_array_samespeciesC(5), &
        pot_samespeciesC(5),sec_der_samespeciesC(5), &
        r_array_thiFeC(4),thiFeC(4),sec_der_thiFeC(4)

    real(8), allocatable:: &
        secder_smspc1(:), &
        secder_smspc2(:),r_potdifspc(:),secder_difspc(:), &
        rho_embFe(:),embFe(:),secder_embFe(:), &
        rho_embC(:),embC(:),secder_embC(:) !For external look-up tables

    logical thiaccptindepndt
    integer embfunctype
    integer, allocatable:: typethi(:,:),typepairpot(:,:)
    !Variables for splines:
    integer narrpairpotXX,narrthiXX
    real(8), allocatable:: r_pairpotXXarr(:),pairpotXXarr(:), &
        secderpairpotXX(:),r_thiXXarr(:),thiXXarr(:), &
        secderthiXX(:)
    !---------------------------
end module m_meamparameters
