module m_neighborlist
    !List of coordinates and neighbors
    integer n_inequivalentsites,nnatoms
    integer, allocatable:: species(:),n_neighbors(:),neighborlist(:,:)
    real(8), allocatable:: xyz(:,:)
end module m_neighborlist
