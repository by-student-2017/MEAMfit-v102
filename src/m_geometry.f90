module m_geometry
    !List of coordinates and neighbors
    integer nstruct,maxnnatoms,istr,smlsepnperspcStr_numEntries
    real(8) rij,rij2,rij3
    integer, allocatable:: gn_inequivalentsites(:), & !struct
        gspecies(:,:),        & !nnatom, struct
        gn_neighbors(:,:),    & !site, struct
        gneighborlist(:,:,:), & !neighborsite,site,struct,
        gn_forces(:),gn_C(:)    !gn_C is the no. of C atoms
    logical fastForce
    real(8), allocatable:: gxyz(:,:,:),diststr(:,:,:,:,:), &
        dist2str(:,:,:,:,:),dist3str(:,:,:,:,:), &
        dxstr(:,:,:,:,:),dystr(:,:,:,:,:), &
        dzstr(:,:,:,:,:), & !xyz,nnatom,struct
        gxyz_backup(:,:,:),smallestsepnStr(:), &
        smlsepnperspcStr(:,:), &
        lrgsepnperspcStr(:,:), &
        smlsepnperspc_overall(:,:), &
        lrgsepnperspc_overall(:,:), &
        sdSepn(:,:), &
        avgSepn(:,:), &
        nSepn(:,:)

contains
    real(8) function distance(i,j)
        integer i,j
        distance=sqrt( (gxyz(1,i,istr)-gxyz(1,j,istr))**2+ &
            (gxyz(2,i,istr)-gxyz(2,j,istr))**2+ &
            (gxyz(3,i,istr)-gxyz(3,j,istr))**2 )
    end function distance
    real(8) function distance2(i,j)
        integer i,j
        distance2=(gxyz(1,i,istr)-gxyz(1,j,istr))**2+ &
            (gxyz(2,i,istr)-gxyz(2,j,istr))**2+ &
            (gxyz(3,i,istr)-gxyz(3,j,istr))**2
    end function distance2
end module m_geometry
