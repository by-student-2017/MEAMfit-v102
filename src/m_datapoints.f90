module m_datapoints
    !Input energy or force data (e.g. from DFT)
    logical, allocatable:: ensureminimum(:),freeenergy(:)
    integer, allocatable:: optimizeforce(:),optimizeforce_backup(:)
    integer, allocatable:: rlxstruc(:)
    integer ndatapoints       !no. of datapoints in file
    real(8) V1,V2,V3,e1,e2,e3,e4,e5,e6,e7,e8,e9, &
        gam4,gam5,gam6,gam7,gam8,gam9, &
        V(21)
    real(8), allocatable:: truedata(:) !size: ndatapoints
    logical, allocatable:: optforce(:,:) !size: maxatoms,nposcars
end module m_datapoints
