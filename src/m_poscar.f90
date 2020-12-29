module m_poscar
    !Data from the poscar file
    integer nspecies,natoms,norigatoms,natoms_old
    integer, allocatable:: z(:),nat(:),zz(:)
    real(8) tr(3,3),trinv(3,3)
    real(8), allocatable:: coordinates(:,:),coords_old(:,:)
end module m_poscar
